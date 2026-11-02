# ══════════════════════════════════════════════════════════════════════════
# EKONOMETRIA PANELOWA — DETERMINANTY OCZEKIWANEJ DŁUGOŚCI ŻYCIA
# Próba: 37 krajów OECD, lata 2010–2019
# ══════════════════════════════════════════════════════════════════════════

library(tidyverse)
library(plm)
library(lmtest)
library(sandwich)
library(car)
library(nortest)
library(ggcorrplot)
library(stargazer)
library(gridExtra)

# ── Wczytanie i przygotowanie danych ──────────────────────────────────────

df_raw <- read.csv2("dane.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)

colnames(df_raw) <- c("country", "year",
                      "life_exp", "gdp", "health_exp",
                      "doctors", "beds", "inf_mort",
                      "obesity", "smoking",
                      "alcohol", "education")

df <- df_raw %>%
  mutate(
    year    = as.integer(year),
    country = as.factor(country),
    across(where(is.character), as.numeric)
  ) %>%
  filter(year <= 2019)

# Transformacje log
df <- df %>%
  mutate(
    log_gdp        = log(gdp),
    log_health_exp = log(health_exp),
    log_inf_mort   = log(inf_mort)
  )

# Centrowanie doctors (globalne i within)
df <- df %>%
  mutate(
    doctors_c  = doctors - mean(doctors),
    obesity_c  = obesity - mean(obesity)
  ) %>%
  group_by(country) %>%
  mutate(doctors_wc = doctors - mean(doctors)) %>%
  ungroup()

# Sprawdzenie struktury
cat("=== STRUKTURA DANYCH ===\n")
cat("Liczba krajów:", n_distinct(df$country), "\n")
cat("Lata:", min(df$year), "-", max(df$year), "\n")
cat("Liczba obserwacji:", nrow(df), "\n")
cat("\nBraki danych:\n")
print(colSums(is.na(df)))

# ══════════════════════════════════════════════════════════════════════════
# KROK 2: STATYSTYKI OPISOWE I WYKRESY WSTĘPNE
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== STATYSTYKI OPISOWE ===\n")
df %>%
  select(life_exp, log_gdp, log_inf_mort, doctors, obesity,
         smoking, alcohol, education) %>%
  summary() %>%
  print()

cat("\n=== DEKOMPOZYCJA WARIANCJI life_exp ===\n")
df_var <- df %>%
  group_by(country) %>%
  mutate(life_exp_mean = mean(life_exp)) %>%
  ungroup()

cat("Wariancja OGÓŁEM:         ", round(var(df$life_exp), 4), "\n")
cat("Wariancja MIĘDZY krajami: ", round(var(df_var$life_exp_mean), 4), "\n")
cat("Wariancja WEWNĄTRZ kraju: ", round(mean(tapply(df$life_exp, df$country, var)), 4), "\n")

dir.create("wykresy", showWarnings = FALSE)

# ── Wykres 1: Trajektorie life_exp w czasie ───────────────────────────────
png("wykresy/01_trajektorie.png", width = 1200, height = 700, res = 120)
ggplot(df, aes(x = year, y = life_exp, group = country)) +
  geom_line(alpha = 0.4, color = "steelblue") +
  geom_smooth(aes(group = 1), method = "loess", color = "red",
              linewidth = 1.2, se = FALSE) +
  theme_minimal() +
  labs(title = "Oczekiwana długość życia w czasie (37 krajów OECD, 2010–2019)",
       x = "Rok", y = "Oczekiwana długość życia (lata)")
dev.off()

# ── Wykres 2: Boxplot według kraju ────────────────────────────────────────
png("wykresy/02_boxplot_kraje.png", width = 1000, height = 900, res = 120)
ggplot(df, aes(x = reorder(country, life_exp, median), y = life_exp)) +
  geom_boxplot(fill = "steelblue", alpha = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Zróżnicowanie oczekiwanej długości życia według kraju",
       x = "", y = "Oczekiwana długość życia (lata)")
dev.off()

# ── Wykres 3: Heterogeniczność między krajami (trajektorie kolorowane) ────
png("wykresy/03_heterogenicznosc.png", width = 1200, height = 700, res = 120)
df %>%
  group_by(country) %>%
  mutate(country_mean = mean(life_exp)) %>%
  ungroup() %>%
  mutate(country = reorder(country, country_mean)) %>%
  ggplot(aes(x = year, y = life_exp, group = country, color = country_mean)) +
  geom_line(linewidth = 0.7) +
  scale_color_gradient(low = "steelblue", high = "tomato", name = "Średnia\nkraju") +
  theme_minimal() +
  labs(title = "Trajektorie life_exp — widoczna heterogeniczność między krajami",
       subtitle = "Kolor = średni poziom kraju; nachylenie linii = zmienność wewnątrzkrajowa",
       x = "Rok", y = "Oczekiwana długość życia (lata)")
dev.off()

# ── Wykres 4: Histogramy — uzasadnienie transformacji log ─────────────────
png("wykresy/04_histogramy_log.png", width = 1200, height = 800, res = 120)
par(mfrow = c(2, 3))
hist(df$gdp,            main = "gdp (surowe)",        xlab = "", col = "steelblue")
hist(df$log_gdp,        main = "log_gdp",              xlab = "", col = "tomato")
hist(df$inf_mort,       main = "inf_mort (surowe)",    xlab = "", col = "steelblue")
hist(df$log_inf_mort,   main = "log_inf_mort",         xlab = "", col = "tomato")
hist(df$health_exp,     main = "health_exp (surowe)",  xlab = "", col = "steelblue")
hist(df$log_health_exp, main = "log_health_exp",       xlab = "", col = "tomato")
par(mfrow = c(1, 1))
dev.off()

# ══════════════════════════════════════════════════════════════════════════
# KROK 3: MACIERZ KORELACJI
# ══════════════════════════════════════════════════════════════════════════

cor_matrix <- df %>%
  select(life_exp, log_gdp, log_inf_mort, doctors,
         obesity, smoking, alcohol, education) %>%
  cor(use = "complete.obs")

cat("\n=== MACIERZ KORELACJI ===\n")
print(round(cor_matrix, 3))
cat("\nKorelacja log_gdp ~ log_health_exp:",
    round(cor(df$log_gdp, df$log_health_exp), 3), "\n")

# ── Wykres 5: Macierz korelacji ───────────────────────────────────────────
png("wykresy/05_macierz_korelacji.png", width = 1000, height = 900, res = 120)
ggcorrplot(cor_matrix,
           method = "square",
           type = "lower",
           lab = TRUE,
           lab_size = 3.5,
           colors = c("#d73027", "white", "#1a9850"),
           title = "Macierz korelacji (zmienne modelu finalnego)")
dev.off()

# ══════════════════════════════════════════════════════════════════════════
# KROK 4a: OLS PEŁNY
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== OLS PEŁNY ===\n")
ols_full <- lm(life_exp ~ log_gdp + log_inf_mort + doctors + beds +
                 obesity + smoking + alcohol + education,
               data = df)
summary(ols_full)

cat("\nVIF:\n")
vif(ols_full)

cat("\nTest normalności reszt (Lilliefors):\n")
lillie.test(residuals(ols_full))

cat("\nTest heteroskedastyczności (Breusch-Pagan):\n")
bptest(ols_full)

cat("\nTest autokorelacji (Durbin-Watson):\n")
dwtest(ols_full)

cat("\nTest RESET:\n")
resettest(ols_full, power = 2:3, type = "fitted")

# ── Wykres 6: Diagnostyka OLS ─────────────────────────────────────────────
png("wykresy/06_diagnostyka_ols.png", width = 1200, height = 1000, res = 120)
par(mfrow = c(2, 2))
plot(ols_full, main = "Diagnostyka OLS")
par(mfrow = c(1, 1))
dev.off()

# ══════════════════════════════════════════════════════════════════════════
# KROK 4b: OLS BEZ BEDS
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== OLS ZREDUKOWANY (bez beds) ===\n")
ols_red <- lm(life_exp ~ log_gdp + log_inf_mort + doctors +
                obesity + smoking + alcohol + education,
              data = df)
summary(ols_red)

cat("\nTest F: czy beds można usunąć?\n")
anova(ols_red, ols_full)

# ══════════════════════════════════════════════════════════════════════════
# KROK 4c: OLS Z KWADRATAMI (naprawa formy funkcyjnej)
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== OLS Z KWADRATAMI (po centrowaniu) ===\n")
ols_quad_c <- lm(life_exp ~ log_gdp + log_inf_mort +
                   doctors_c + I(doctors_c^2) +
                   obesity_c + I(obesity_c^2) +
                   smoking + alcohol + education,
                 data = df)
summary(ols_quad_c)

cat("\nRESET po centrowaniu:\n")
resettest(ols_quad_c, power = 2:3, type = "fitted")

cat("\nVIF po centrowaniu:\n")
vif(ols_quad_c)

cat("\nTest F: ols_quad_c vs ols_red:\n")
anova(ols_red, ols_quad_c)

# ══════════════════════════════════════════════════════════════════════════
# KROK 5: MODELE PANELOWE
# ══════════════════════════════════════════════════════════════════════════

pdf_data <- pdata.frame(df, index = c("country", "year"))

# FE liniowy
cat("\n=== MODEL FE LINIOWY ===\n")
fe_full <- plm(life_exp ~ log_gdp + log_inf_mort + doctors +
                 obesity + smoking + alcohol + education,
               data = pdf_data,
               model = "within")
summary(fe_full)

# FE z nieliniowością doctors
cat("\n=== MODEL FE Z NIELINIOWOŚCIĄ DOCTORS ===\n")
fe_quad <- plm(life_exp ~ log_gdp + log_inf_mort +
                 doctors_wc + I(doctors_wc^2) +
                 obesity + smoking + alcohol + education,
               data = pdf_data,
               model = "within")
summary(fe_quad)

cat("\nTest F: czy doctors² istotny w FE?\n")
pFtest(fe_quad, fe_full)

# RE
cat("\n=== MODEL RE ===\n")
re_quad <- plm(life_exp ~ log_gdp + log_inf_mort +
                 doctors_wc + I(doctors_wc^2) +
                 obesity + smoking + alcohol + education,
               data = pdf_data,
               model = "random")
summary(re_quad)

# ══════════════════════════════════════════════════════════════════════════
# KROK 6: TESTY WYBORU MODELU
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== TEST BREUSCHA-PAGANA (RE vs OLS) ===\n")
plmtest(re_quad, type = "bp")

cat("\n=== TEST HAUSMANA (FE vs RE) ===\n")
phtest(fe_quad, re_quad)

# ══════════════════════════════════════════════════════════════════════════
# KROK 7: DIAGNOSTYKA MODELU FINALNEGO FE
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== DIAGNOSTYKA FE FINALNEGO ===\n")

cat("\nAutokorelacja (Wooldridge):\n")
pbgtest(fe_quad)

cat("\nHeteroskedastyczność (Breusch-Pagan):\n")
bptest(fe_quad)

cat("\nZależność przekrojowa (Pesaran):\n")
pcdtest(fe_quad)

# ══════════════════════════════════════════════════════════════════════════
# KROK 8: BŁĘDY DRISCOLL-KRAAY + PUNKT MAKSIMUM
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== FE FINALNY Z BŁĘDAMI DRISCOLL-KRAAY ===\n")
coeftest(fe_quad, vcov = vcovSCC(fe_quad, type = "HC3", cluster = "time"))

# Punkt maksimum dla doctors
max_doctors_dev   <- coef(fe_quad)["doctors_wc"] / (2 * abs(coef(fe_quad)["I(doctors_wc^2)"]))
max_doctors_level <- mean(df$doctors) + max_doctors_dev
cat("\nPunkt maksimum doctors:", round(max_doctors_level, 3), "lekarzy/1000\n")
cat("Odsetek obserwacji powyżej maksimum:",
    round(mean(df$doctors > max_doctors_level) * 100, 1), "%\n")

# Tabela: które kraje są po której stronie krzywej
cat("\n=== KRAJE VS PUNKT MAKSIMUM ===\n")
df %>%
  group_by(country) %>%
  summarise(avg_doctors = round(mean(doctors), 2)) %>%
  mutate(
    vs_max = round(avg_doctors - max_doctors_level, 2),
    status = ifelse(avg_doctors < max_doctors_level, "poniżej maksimum", "powyżej maksimum")
  ) %>%
  arrange(avg_doctors) %>%
  print(n = 37)

# ── Wykres 7: Odwrócona U dla doctors ────────────────────────────────────
png("wykresy/07_odwrocona_u_doctors.png", width = 900, height = 700, res = 120)

b1 <- coef(fe_quad)["doctors_wc"]
b2 <- coef(fe_quad)["I(doctors_wc^2)"]
doctors_seq <- seq(min(df$doctors), max(df$doctors), length.out = 200)
doctors_wc_seq <- doctors_seq - mean(df$doctors)
y_seq <- b1 * doctors_wc_seq + b2 * doctors_wc_seq^2

curve_df <- data.frame(doctors = doctors_seq, effect = y_seq)

country_means <- df %>%
  group_by(country) %>%
  summarise(avg_doctors = mean(doctors))

ggplot(curve_df, aes(x = doctors, y = effect)) +
  geom_line(color = "steelblue", linewidth = 1.3) +
  geom_vline(xintercept = max_doctors_level, color = "red",
             linetype = "dashed", linewidth = 1) +
  annotate("text", x = max_doctors_level + 0.1, y = max(y_seq) * 0.95,
           label = paste0("Maksimum\n", round(max_doctors_level, 2), " lek./1000"),
           color = "red", hjust = 0, size = 3.5) +
  geom_rug(data = country_means, aes(x = avg_doctors, y = NULL),
           color = "gray40", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Nieliniowy efekt liczby lekarzy na oczekiwaną długość życia",
       subtitle = "Odwrócona U — efekt FE (within), błędy Driscoll-Kraay",
       x = "Liczba lekarzy na 1000 mieszkańców",
       y = "Efekt cząstkowy na life_exp (lata)")
dev.off()

# ── Wykres 8: Efekty stałe krajów (fixed effects) ────────────────────────
png("wykresy/08_efekty_stale.png", width = 1000, height = 900, res = 120)
fixef(fe_quad) %>%
  as.data.frame() %>%
  rownames_to_column("country") %>%
  rename(effect = ".") %>%
  arrange(effect) %>%
  mutate(country = factor(country, levels = country)) %>%
  ggplot(aes(x = effect, y = country)) +
  geom_col(aes(fill = effect > 0), show.legend = FALSE) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Efekty stałe krajów (fixed effects)",
       subtitle = "Nieobserwowalna heterogeniczność między krajami",
       x = "Efekt stały (odchylenie od średniej)", y = "")
dev.off()

# ══════════════════════════════════════════════════════════════════════════
# KROK 9: TABELA PORÓWNAWCZA OLS / RE / FE
# ══════════════════════════════════════════════════════════════════════════

se_ols   <- sqrt(diag(vcovHC(ols_full, type = "HC3")))
se_re_dk <- sqrt(diag(vcovSCC(re_quad, type = "HC3", cluster = "time")))
se_fe_dk <- sqrt(diag(vcovSCC(fe_quad, type = "HC3", cluster = "time")))

stargazer(ols_full, re_quad, fe_quad,
          type = "text",
          column.labels = c("OLS (HC3)", "RE (Driscoll-Kraay)", "FE (Driscoll-Kraay)"),
          title = "Porównanie modeli — OLS, RE, FE",
          digits = 3,
          se = list(se_ols, se_re_dk, se_fe_dk),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = c("OLS: błędy HC3; RE i FE: błędy Driscoll-Kraay (1998)",
                    "Test Hausmana: FE preferowany nad RE (p = 0.0008)"))

# ══════════════════════════════════════════════════════════════════════════
# KROK 10: ROBUSTNESS CHECK — TWO-WAY FE
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== ROBUSTNESS CHECK: TWO-WAY FE ===\n")
fe_twoway <- plm(life_exp ~ log_gdp + log_inf_mort +
                   doctors_wc + I(doctors_wc^2) +
                   obesity + smoking + alcohol + education,
                 data = pdf_data,
                 model = "within",
                 effect = "twoways")

coeftest(fe_twoway, vcov = vcovSCC(fe_twoway, type = "HC3", cluster = "time"))

se_fe_dk     <- sqrt(diag(vcovSCC(fe_quad,   type = "HC3", cluster = "time")))
se_twoway_dk <- sqrt(diag(vcovSCC(fe_twoway, type = "HC3", cluster = "time")))

stargazer(fe_quad, fe_twoway,
          type = "text",
          column.labels = c("FE one-way (główny)", "FE two-way (robustness)"),
          title = "Robustness check — Two-way FE",
          digits = 3,
          se = list(se_fe_dk, se_twoway_dk),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = "Błędy Driscoll-Kraay (1998)")

# ══════════════════════════════════════════════════════════════════════════
# KROK 11: DODATKOWE WYKRESY DIAGNOSTYCZNE I INTERPRETACYJNE
# ══════════════════════════════════════════════════════════════════════════

# ── Wykres 9: Scatter — PKB vs life_exp (krzywa Prestona) ────────────────
png("wykresy/09_preston_curve.png", width = 900, height = 700, res = 120)
df %>%
  group_by(country) %>%
  summarise(avg_gdp = mean(gdp), avg_life_exp = mean(life_exp)) %>%
  ggplot(aes(x = avg_gdp, y = avg_life_exp, label = country)) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(vjust = -0.7, size = 2.5) +
  geom_smooth(method = "loess", color = "red", se = TRUE, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Krzywa Prestona — PKB per capita vs oczekiwana długość życia",
       subtitle = "Kraje OECD, średnie 2010–2019",
       x = "PKB per capita (USD)", y = "Oczekiwana długość życia (lata)")
dev.off()

# ── Wykres 10: Scatter — education vs life_exp ───────────────────────────
png("wykresy/10_education_vs_lifeexp.png", width = 900, height = 700, res = 120)
df %>%
  group_by(country) %>%
  summarise(avg_edu = mean(education), avg_life_exp = mean(life_exp)) %>%
  ggplot(aes(x = avg_edu, y = avg_life_exp, label = country)) +
  geom_point(color = "steelblue", size = 3) +
  geom_text(vjust = -0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "red", se = TRUE, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Edukacja vs oczekiwana długość życia",
       subtitle = "Kraje OECD, średnie 2010–2019",
       x = "Średnia liczba lat edukacji", y = "Oczekiwana długość życia (lata)")
dev.off()

# ── Wykres 11: Scatter — inf_mort vs life_exp ────────────────────────────
png("wykresy/11_infmort_vs_lifeexp.png", width = 900, height = 700, res = 120)
df %>%
  group_by(country) %>%
  summarise(avg_inf = mean(inf_mort), avg_life_exp = mean(life_exp)) %>%
  ggplot(aes(x = avg_inf, y = avg_life_exp, label = country)) +
  geom_point(color = "tomato", size = 3) +
  geom_text(vjust = -0.7, size = 2.5) +
  geom_smooth(method = "lm", color = "darkred", se = TRUE, linetype = "dashed") +
  theme_minimal() +
  labs(title = "Umieralność niemowląt vs oczekiwana długość życia",
       subtitle = "Kraje OECD, średnie 2010–2019",
       x = "Umieralność niemowląt (na 1000 urodzeń)", y = "Oczekiwana długość życia (lata)")
dev.off()

# ── Wykres 12: Reszty FE vs dopasowane wartości ───────────────────────────
png("wykresy/12_reszty_fe.png", width = 900, height = 700, res = 120)
data.frame(
  fitted = fitted(fe_quad),
  resid  = residuals(fe_quad)
) %>%
  ggplot(aes(x = fitted, y = resid)) +
  geom_point(alpha = 0.4, color = "steelblue") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", color = "orange", se = FALSE, linewidth = 0.8) +
  theme_minimal() +
  labs(title = "Reszty modelu FE vs dopasowane wartości",
       x = "Wartości dopasowane", y = "Reszty")
dev.off()

cat("\n=== WSZYSTKIE WYKRESY ZAPISANE W FOLDERZE wykresy/ ===\n")
cat("01_trajektorie.png\n")
cat("02_boxplot_kraje.png\n")
cat("03_heterogenicznosc.png\n")
cat("04_histogramy_log.png\n")
cat("05_macierz_korelacji.png\n")
cat("06_diagnostyka_ols.png\n")
cat("07_odwrocona_u_doctors.png\n")
cat("08_efekty_stale.png\n")
cat("09_preston_curve.png\n")
cat("10_education_vs_lifeexp.png\n")
cat("11_infmort_vs_lifeexp.png\n")
cat("12_reszty_fe.png\n")

png("wykresy/03_heterogenicznosc.png", width = 1200, height = 700, res = 120)
df %>%
  group_by(country) %>%
  mutate(country_mean = mean(life_exp)) %>%
  ungroup() %>%
  mutate(country = reorder(country, country_mean)) %>%
  ggplot(aes(x = year, y = life_exp, group = country, color = country_mean)) +
  geom_line(linewidth = 0.7) +
  scale_color_gradient(low = "steelblue", high = "tomato", name = "Średnia\nkraju") +
  scale_x_continuous(breaks = 2010:2019) +
  theme_minimal() +
  labs(title = "Trajektorie oczekiwanej długości życia (37 krajów OECD, 2010–2019)",
       subtitle = "Kolor = średni poziom kraju",
       x = "Rok", y = "Oczekiwana długość życia (lata)")
dev.off()

se_fe_dk     <- sqrt(diag(vcovSCC(fe_quad,   type = "HC3", cluster = "time")))
se_twoway_dk <- sqrt(diag(vcovSCC(fe_twoway, type = "HC3", cluster = "time")))

stargazer(fe_quad, fe_twoway,
          type = "text",
          column.labels = c("FE one-way (główny)", "FE two-way (robustness)"),
          title = "Tabela 2. Sprawdzian odporności",
          digits = 3,
          se = list(se_fe_dk, se_twoway_dk),
          star.cutoffs = c(0.05, 0.01, 0.001),
          notes = "Błędy Driscoll-Kraay (1998) w nawiasach.")
