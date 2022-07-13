
# MODEL VAR, cauzalitate Granger ------------------------------------------

# Loading the necessary libraries -----------------------------------------
library(urca)
library(vars)
library(mFilter)
library(tseries)
library(forecast)
library(tidyverse)
library(stargazer)
library(seasonal)

# BTC return vs Google Trends ---------------------------------------------

# Loading the data --------------------------------------------------------
bd <- read.csv(file.choose())
bd = bd[-c(1),]
View(bd)

# Graphs ------------------------------------------------------------------

# Scatterplot
ggplot(data = bd) + 
  geom_point(mapping = aes(x = BTC_Close , y = est_hits)) +
  xlab('Pretul de inchidere Bitcoin ') +
  ylab('Google Trends') + 
  ggtitle('Norul de puncte dintre Google Trends si randamentul pretului de inchidere al BTC')+
  theme_bw()

# Declare variables as ts 
BTC_return <- ts(bd$Randament_BTC, start = c(2011,9,2),frequency = 365)
BTC_Close = ts(bd$BTC_Close, start = c(2011,9,2), frequency = 365)
google <-ts(bd$est_hits, start = c(2011,9,2),frequency = 365)

# Series graphs
autoplot(cbind(log(google), BTC_return)) +
  ylab('') +
  ggtitle('Graficul seriilor') +
  theme_bw()

# Avand in vedere ca datele au ordine de marime diferite, pentru a vizualiza grafic seriile, se vor standardiza. 
# FCT: scale(x, center=TRUE, scale=TRUE)
BTC_Close_scaled = scale(BTC_Close, center=TRUE, scale=TRUE)
BTC_return_scaled = scale(BTC_return, center=TRUE, scale=TRUE)
google_scaled = scale(google, center=TRUE, scale=TRUE)

google_scaled = as.ts(google_scaled)
BTC_return_scaled = as.ts(BTC_return_scaled)
BTC_Close_scaled = as.ts(BTC_Close_scaled)

autoplot(cbind(BTC_Close_scaled, google_scaled)) +
  ylab('') +
  ggtitle('Graficul seriilor - standardizate ') +
  theme_bw()

# Deseasonalization -> Seasonal Adjustment of Daily Time Series -----------

# -------STL DECOMPOSITION
#--------------------BITCOIN CLOSE PRICE
BTC_Close_decompose <- stl(BTC_Close, s.window = 5) 
autoplot(BTC_Close_decompose) +
  ggtitle("STL decomposition of daily Bitcoin close price")

#remove sesonality
BTC_Close_season_adj <- seasadj(BTC_Close_decompose)
plot(BTC_Close_season_adj)

autoplot(BTC_Close, series="Data") +
  autolayer(trendcycle(BTC_Close_decompose), series="Trend") +
  autolayer(seasadj(BTC_Close_decompose), series="Seasonally Adjusted") +
  xlab("Year") + ylab("$") +
  ggtitle("Daily Bitcoin Close Price $") +
  scale_colour_manual(values=c("black","blue","red"),
                      breaks=c("Data","Seasonally Adjusted","Trend"))  
Box.test(abs(BTC_Close_season_adj),lag = 12, type='Ljung')

# ----BITCOIN RETURN

BTC_return_decompose <- stl(BTC_return, s.window = 5) 
autoplot(BTC_return_decompose) +
  ggtitle("STL decomposition of daily Bitcoin return")

#remove sesonality
BTC_return_season_adj <- seasadj(BTC_return_decompose)
plot(BTC_return_season_adj)

autoplot(BTC_return, series="Data") +
  autolayer(trendcycle(BTC_return_decompose), series="Trend") +
  autolayer(seasadj(BTC_return_decompose), series="Seasonally Adjusted") +
  xlab("Year") + ylab("%") +
  ggtitle("daily bitcoin return-%") +
  scale_colour_manual(values=c("black","blue","red"),
                      breaks=c("Data","Seasonally Adjusted","Trend"))  
Box.test(abs(BTC_return_season_adj),lag = 12, type='Ljung')

# ----GOOGLE TRENDS

google_decompose <- stl(google, s.window = 5) 
autoplot(google_decompose) +
  ggtitle("STL decomposition of daily no of Google searches for Bitcoin - Google Trends")

#remove sesonality
google_season_adj <- seasadj(google_decompose)
plot(google_season_adj)

autoplot(google, series="Data") +
  autolayer(trendcycle(google_decompose), series="Trend") +
  autolayer(seasadj(google_decompose), series="Seasonally Adjusted") +
  xlab("Year") + ylab("%") +
  ggtitle("daily no of Google searches for Bitcoin - Google Trends") +
  scale_colour_manual(values=c("black","blue","red"),
                      breaks=c("Data","Seasonally Adjusted","Trend"))  
Box.test(abs(google_season_adj),lag = 12, type='Ljung')



# Q: Series is stationary or non-stationary? ------------------------------
BTC_return = BTC_return_season_adj 
google = google_season_adj

# Determinarea persistentei modelului
ggtsdisplay(BTC_return) 
ggtsdisplay(google)

# Unit-Root Tests (we choose the most complex variant of the ADF Test)
# ---------GOOGLE
# -----ADF
adf.google <- ur.df(google, type = "trend", selectlags = "AIC")
summary(adf.google) # serie stationara

adf.google <- ur.df(google, type = "drift", selectlags = "AIC")
summary(adf.google) # serie stationara

adf.google <- ur.df(google, type = "none", selectlags = "AIC")
summary(adf.google) # serie stationara

#-----KPSS
kpss_google <- google %>%
  ur.kpss(., type="mu", use.lag=NULL) # mu for level stationary
summary(kpss_google)     #non-stationary   

#-----ZA TEST
za_google <- google %>%
  ur.za(., model="trend", lag = 1)
summary(za_google)
plot(za_google) 

# ---------BITCOIN RETURN
adf.BTC_return <- ur.df(BTC_return, type = "trend", selectlags = "AIC")
summary(adf.BTC_return) # serie stationara 

adf.BTC_return <- ur.df(BTC_return, type = "drift", selectlags = "AIC")
summary(adf.BTC_return) # serie stationara 

adf.BTC_return <- ur.df(BTC_return, type = "none", selectlags = "AIC")
summary(adf.BTC_return) # serie stationara 

#-----KPSS
kpss_btc <- BTC_return %>%
  ur.kpss(., type="mu", use.lag=NULL) # mu for level stationary
summary(kpss_btc)     #stationary, accept the null hypothesis, t stat < all the critical values => series is I(0)

#-----ZA TEST
za_btc <- BTC_return %>%
  ur.za(., model="trend", lag = 1)
summary(za_btc)
plot(za_btc) 

# VAR Methodology ---------------------------------------------------------

# Finding the optimal lag -------------------------------------------------
df <- cbind(BTC_return,google)
colnames(df) <- cbind('BTC_return','Google_Trends')
View(df)
lagselect <- VARselect(df,lag.max = 14, type = 'const')
lagselect
lagselect$selection # lagul 7, conform HQ, lagul 6, conform SC, lagul 24, conform AIC si FPE 


# VAR ---------------------------------------------------------------------
model1 <- VAR(df, p = 6, type = 'const', season = NULL, exog = NULL)
summary(model1)

model2 <- VAR(df, p = 7, type = 'const', season = NULL, exog = NULL)
summary(model2)

model3 <- VAR(df, p = 14, type = 'const', season = NULL, exog = NULL)
summary(model3)

# Radacinile unitate (roots of the characteristic polynomial) < 1
# toate se afla in interiorul cercului unitate => model stabil
# Ecuatia BTC
# Ecuatia Google Trends 
# Reziduurile nu se coreleaza BTC ~ Google Trends (0.02436 corelatie slaba)

# O alta modalitate de afisare a rezultatelor modelului 
stargazer(model2[['varresult']], type = 'text')


# Residuals diagnostics ---------------------------------------------------

# H1: Autocorrelation -----------------------------------------------------

Serial2 <- serial.test(model2, lags.pt = 9, type = 'PT.asymptotic')
Serial2 # pvalue < 0.1  avem autocorelare in reziduuri

# H2: Homoskedasticity ----------------------------------------------------

Arch2 <- vars::arch.test(model2,lags.multi = 9,multivariate.only = TRUE)
Arch2 # pvalue < 0.05 modelul prezinta heteroschedasticitate la 95%

# H3: Normality -----------------------------------------------------------

Norm2 <- normality.test(model2, multivariate.only = TRUE)
Norm2 # pvalue JB < 0.05 reziduurile nu sunt normal distribuite


# Structural breaks in series ---------------------------------------------
# Testarea pentru rupturi in serie

Stability2 <- stability(model2,type = 'OLS-CUSUM')
plot(Stability2) # model aproape stabil, seria BTC nu depaseste intervalul rosu, dar seria Google depaseste putin intervalul rosu

# Granger Causality -------------------------------------------------------

# Cauzalitatea Granger determina daca un model care utilizeaza valorile trecute
# si prezente ale lui X si valorile trecute si prezente ale lui Y prezinta erori
# mai mici de prognozare fata de un model care nu prezinta cauzalitate Granger
# Cu alte cuvinte cauzalitatea Granger raspunde la intrebarea: trecutul variabilei X
# ajuta la imbunatatirea predictiei valorilor lui y?

# Pentru a testa cauzalitatea Granger, trebuie sa ne asiguram ca seriile noastre sunt
# stationare si ca nu avem autocorelare

# Ipotezele cauzalitatii Granger
# H0: valorile cu lag ale lui X, nu explica variatia in Y 
# H1: valorile cu lag ale lui X, explica variatia in Y

# H0: variabila X nu prezinta cauzalitate Granger pentru variabila Y
# H1: Variabila X prezinta cauzalitate Granger pentru variabila Y

# Cauzalitatea Granger construita in jurul testului F isi poate pierde din putere
# atunci cand avem un numar mare de variabile si laguri
# Varianta alternativa a cauzalitatii Granger folosita in cazul in care avem multe laguri
# si un numar mare de serii in modelul VAR si se bazeaza pe testul Wald si distributia Chi-patrat
# Cu alte cuvinte cauzalitatea Granger-Wald denumita si cauzalitate instantanee
# raspunde la intrebarea: cunoasterea viitorului lui x ma ajuta sa prezic mai bine 
# viitorul lui y?

# Ipotezele pentru ambele teste de cauzalitate Granger raman la fel, in output 
# cea clasica o gasim sub forma de Granger, iar cea bazata pe Wald o gasim 
# sub forma Instant

GrangerBTC2 <- causality(model2, cause = 'BTC_return')
GrangerBTC2 # p < 0.1, p-value = 0.00411 => BTC prezinta cauzalitate Granger cu Google Trends
GrangerGoogle2 <- causality(model2, cause = 'Google_Trends')
GrangerGoogle2 # p < 0.1, p-value = 0.03009 => Google Trends prezinta cauzalitate Granger cu BTC

# Impulse-Response Function (IRF) -----------------------------------------
#-----M2
BTCirf2 <- irf(model2, impulse = 'Google_Trends', response = 'BTC_return', 
               n.ahead = 20, boot = TRUE, ci=0.90) # n.ahead = perioade in viitor
                                          # boot = TRUE pentru a crea intervale de incredere
                                          # ci - nivelul de semnificatie
plot(BTCirf2, ylab = 'BTC_return', main = 'Răspunsul BTC la șocurile Google Trends')

Googleirf2 <- irf(model2, impulse = 'BTC_return', response = 'Google_Trends', 
                 n.ahead = 20, boot = TRUE, ci=0.90)
plot(Googleirf2, ylab = 'Google_Trends', main = 'Răspunsul Google Trends la șocurile BTC')

# Variance Decomposition --------------------------------------------------

FEVD <- fevd(model2, n.ahead = 20)
plot(FEVD) # graficul ne spune procentul de unde vine socul variabilei
# Pentru ambele variabile socul vin de la variabila in sine si mai putin de la cealalta variabila din model 
FEVD 

# Forecasting -> VAR ------------------------------------------------------

forecast2 <- predict(model2, n.ahead = 100, ci = 0.90) # prognoza pe 100 zile 
forecast2
plot(forecast2, name = 'BTC_return')
plot(forecast2, name = 'Google_Trends')
fanchart(forecast2, names='BTC_return')
fanchart(forecast2, names='Google_Trends')



