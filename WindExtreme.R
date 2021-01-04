## ----setup, echo=FALSE----------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(knitr.purl.inline = TRUE)
knitr::opts_chunk$set(fig.path = "README_figs/README-")


## ----packages, warning=FALSE, message=FALSE-------------------------------------------------------------------------------------
library(readr)
library(dplyr)
library(lubridate) #date
library(xts) #zoo
library(reshape2) # table
library(ggplot2)
require(patchwork) #2 ggplots


## ----read data, warning=FALSE, message=FALSE------------------------------------------------------------------------------------
X545_Porto_FFMax = read_delim("545_Porto_FFMax.csv", 
    ";", escape_double = FALSE, trim_ws = TRUE)

colnames(X545_Porto_FFMax)



## ----set dos parametros de análise----------------------------------------------------------------------------------------------
par(mfrow=c(1,2))

windspeed = X545_Porto_FFMax

u = seq(0,90,1) # extremo do nível u dos excedentes é 90Kms/h

x = numeric()
y = numeric()
en_u = data.frame()

# Critérios de selecção do Ano e demais atributos
ws = subset(windspeed, FF_MAX != -990)
ws = subset(ws, !is.na(FF_MAX))
ws$FF_KM = ws$FF_MAX*3.6# criação de uma nova coluna (FF_KM), com o valor do vento calculado em Kms/h

hist(ws$FF_KM, main="Histograma Ventos",xlab="Kms/h", ylab="Frequency")
boxplot(ws$FF_KM, main="Boxplot Ventos", ylab="Kms/h")
paste("min: ",min(ws$FF_KM),"; max: ",max(ws$FF_KM),"; median: ",median(ws$FF_KM))


## ----máximo por blocos de 14 dias-----------------------------------------------------------------------------------------------

Date.time <- ws %>% select(ANO, MS, DI, HR, MN) %>% mutate(Date_Time = make_datetime(ANO, MS, DI, HR, MN))
ws$Date <- Date.time$Date_Time

ws_2col <- ws[,c(9,8)]; ts.dat <- read.zoo(file = ws_2col); ep <- endpoints(ts.dat,"days", k=14)

FF_KM <- period.apply(x = ts.dat,ep,FUN = max); ws_d14 <- fortify.zoo(FF_KM, name="Date")


## ----representação gráfica------------------------------------------------------------------------------------------------------
plot(ws$Date, ws$FF_KM, xlab="Date", ylab="Wind Velocity (Km/h)",
     cex=1.25, cex.lab=1.25,
     col = "gray", bg = "lightblue", pch=21)
title(main = "Block Maxima Approach by 14 days")
points(ws_d14$Date, ws_d14$FF_KM, col="red", cex=1.5)

m = 157
y_i_sorted = sort(ws_d14$FF_KM)



## ----GEV Maximization-----------------------------------------------------------------------------------------------------------
GEV_M = function(xdat){
  z = list()
  
  #só teremos uma mu, um stdev e um shape
  npmu = 1
  npstdev = 1
  npsh = 1 #shape
  
  mumat = as.matrix(rep(1, length(xdat)))
  sigmat = as.matrix(rep(1, length(xdat)))
  shmat = as.matrix(rep(1, length(xdat)))
  
  #initial values for minimization routine (init for each parameter)
  siginit = sqrt(6 * var(xdat))/pi
  muinit = mean(xdat) - 0.57722 * siginit
  shinit = 0.1
  
  init = c(muinit, siginit, shinit)
  
  gev.MLE = function(a) {
  # computes -log lik of gev model
  mu = identity(mumat %*% a[1])
  sc = identity(sigmat %*% a[2])
	xi = identity(shmat %*% a[3])
	
	y = (xdat - mu)/sc
	y = 1 + xi * y
	#has to be positive
	if(any(y <= 0) || any(sc <= 0)) return(10^6)
	#function 
	sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y) * (1/xi + 1))
  }
  
  #Initial values for the parameters to be optimized over.
  #A function to be minimized (or maximized), with first argument the vector of parameters over which minimization is to take place.
	x = optim(init, gev.MLE, hessian = TRUE, method = "Nelder-Mead",
                   control = list(maxit = 10000))
  #mle
	z$mle = x$par
	
	z$data = xdat

	invisible(z)
}


## ----GEV ML compute-------------------------------------------------------------------------------------------------------------
par_GEV = GEV_M(y_i_sorted)
par_GEV


## ----plot par_GEV---------------------------------------------------------------------------------------------------------------
par(mfrow=c(1,2))

gev.his = function(a, dat){
# Plots histogram of data and fitted density

  h = hist(dat, plot = FALSE)
	x_axis = seq(max(min(h$breaks), (a[1] - a[2]/a[3] + 0.001)), max(h$
			breaks), length = 100)

	# density
	gev.dens = function(a, z){ # evaluates gev density with parameters a at z 
  if(round(a[3],1) != 0)
    (exp( - (1 + (a[3] * (z - a[1]))/a[2])^(-1/a[3])) * (1 + (
			a[3] * (z - a[1]))/a[2])^(-1/a[3] - 1))/a[2]
    else { #gumbel 
      y = (z - a[1])/a[2]
	    (exp( - y) * exp( - exp( - y)))/a[2]
      }
}
	y = gev.dens(a, x_axis)
	#por causa do parâmetro, 
	dGumbel<-function(x,mu,sigma){exp(-((x-mu)/sigma+exp(-(x-mu)/sigma)))/sigma}
	# hist
	hist(dat, freq = FALSE, ylim=c(0,0.04),xlab = "Km/h", ylab = "f(Km/h)", 
		main = "Density Plot")
	# data
	points(dat, rep(0, length(dat)))
	#density
	lines(x_axis,y, col='red')
	#Gumbel
	curve(dGumbel(x, a[1], a[2]), add = TRUE, lwd=2, col="red")
}

gev.his(par_GEV$mle, par_GEV$data)


# gev com os parametros estimados
GEV_func = function(a, z){# gev dist fnc
	exp( - (1 + (a[3] * (z - a[1]))/a[2])^(-1/a[3]))}

GEVf = GEV_func(par_GEV$mle, y_i_sorted[1:(m-1)])
plot(y_i_sorted[1:(m-1)], GEVf)



## ----plot da amostra e dos valores acima de u, para u igual 70------------------------------------------------------------------
ggplot(data=ws, aes(x=Date, y=FF_KM)) + geom_point(aes(colour = cut(FF_KM, c(-Inf, 70, Inf)))) + 
  scale_color_manual(name = "FF_KM",
                     values = c("(-Inf,70]" = "gray",
                                "(70, Inf]" = "red"),
                     labels = c("<70", ">70")) +
  ggtitle("Peak Over Threshold para u = 70 Km/H")


## ----média dos excessos da amostra----------------------------------------------------------------------------------------------
# Função da média dos excessos da amostra - pág. 5 do artigo "The Peak over threshold method for estimating high quantiles of loss distributions"

for (j in 1:length(u)){
  x <- 0
  y <- 0
  for (i in 1: nrow(ws)){
    if(ws$FF_KM[i]>u[j]){
      x <- x + (ws$FF_KM[i]-u[j])
      y <- y + 1
    }
  }
  en_u[j,1] <- u[j]# definição do nível u[j], com registo para os diferentes valores que toma - [0,10,20,..., 90]
  en_u[j,2] <- x # Somatório dos ventos acima de u[j]
  en_u[j,3] <- y # função identidade para os valores dos ventos acima de u[j]
  en_u[j,4] <- x/y # função da média amostral dos excessos
  en_u[j,5] <- (1-y/nrow(ws)) # Quantil... por aqui se ê que temos extremos interessantes e significativos
}


## ----Plot da função da média dos excessos da amostra----------------------------------------------------------------------------
colnames(en_u) <- c("u","x","y","x_y", "quantile")
head(en_u)

plot(en_u$u[1:nrow(en_u)],en_u$x_y, type="l", main="Mean Excess Plot", xlab="Threshold (u)",ylab="Mean Excess") # gráfico da função tal como sugerido pelo gráfico da função da média dos excessos

en_u_cor_10 <- cor(en_u[1:10,4],en_u[1:10,1]) # Cálculo da correlação entre os dados do u e os dados da ME
cor_10 <- lm(en_u[1:10,4]~en_u[1:10,1])
en_u_cor_70 <- cor(en_u[11:70,4],en_u[11:70,1])
cor_70 <- lm(en_u[11:63,4]~en_u[11:63,1])
abline(cor_10,col="red", lty=4)
abline(cor_70,col="blue", lty=2)


## ----Excessos acima de 53Kms/h--------------------------------------------------------------------------------------------------

u_threshold <- 53
ws_u <- subset(ws, ws$FF_KM>u_threshold)
nrow(ws_u)
ws$Date <- as.yearmon(paste(ws$ANO, ws$MS, ws$DI), "%Y %m %d")
ggplot(data=ws, aes(x=Date, y=FF_KM)) + geom_point(aes(colour = cut(FF_KM, c(-Inf, 53, Inf)))) + 
  scale_color_manual(name = "FF_KM",
                     values = c("(-Inf,53]" = "gray",
                                "(53, Inf]" = "blue"),
                     labels = c("<=53", ">53")) +
  ggtitle("Peak Over Threshold para u = 53Km/H")

paste("n.º registos com ventos > 53Kms/h =",nrow(ws_u))


## ----pdf e ecdf dos ventos forte ou extremos------------------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(density(x=(ws_u$FF_KM-u_threshold)), main="pdf de ventos > 53Kms/h", xlab="Ventos Kms/h",ylab="density") # Plot da função densidade para os excessos acima de u (>53Kms/h) 
plot(ecdf(x=(ws_u$FF_KM-u_threshold)), main="ecdf de ventos > 53Kms/h", xlab="Ventos Kms/h",ylab="probability") # f.d. dos excessos acima de u (=70Kms/h)


## ----GPD, warning=FALSE---------------------------------------------------------------------------------------------------------
GPD_M = function(xdat, threshold){
  z = list()
  
  n = length(xdat)
  u = rep(threshold, length.out = n)
  
  xdatu = xdat[xdat > u]
  xind = (1:n)[xdat > u]
	u = u[xind]
	

	#só teremos uma mu, um stdev e um shape
  npmu = 1
  npstdev = 1
  npsh = 1 #shape
  
  mumat = as.matrix(rep(1, length(xdatu)))
  sigmat = as.matrix(rep(1, length(xdatu)))
  shmat = as.matrix(rep(1, length(xdatu)))
  
  
  #initial values for minimization routine (init for each parameter)
  siginit = sqrt(6 * var(xdatu))/pi
  muinit = mean(xdatu, na.rm = TRUE) - 0.57722 * siginit
  shinit = 0.1
  
  init = c(siginit, shinit)
  

  gpd.lik = function(a) {
  # calculates gpd neg log lik
	sc = identity(sigmat %*% (a[seq(1, length = 1)]))
	xi = identity(shmat %*% (a[seq(2, length = 1)]))
	
	y = (xdatu - u)/sc
	y = 1 + xi * y
  l = sum(log(sc)) + sum(log(y) * (1/xi + 1))}

  #Initial values for the parameters to be optimized over.
  #A function to be minimized (or maximized), with first argument the vector of parameters over which minimization is to take place.
  x = optim(init, gpd.lik, hessian = TRUE, method = "Nelder-Mead",
                   control = list(maxit = 10000))
  
  #mle
  z$mu = muinit
	z$mle = x$par
	z$threshold = threshold
	z$nexc = length(xdatu)

	invisible(z)
}

par_GPD = GPD_M(ws$FF_KM, 53)
par_GPD



## ----compute blocks per year----------------------------------------------------------------------------------------------------
#compute blocks per year
block_days = 14
blocks_per_year = 365/block_days; blocks_per_year


## ----function max---------------------------------------------------------------------------------------------------------------
#function
zp = function(mu, stdev, sh, year){#sh - shape
  
  m_blocks = round(blocks_per_year)*year #obter blocos pela transformação
  
  #compute p
  p = 1/m_blocks
  
  yp = -log(1-p)
  if(sh==0){
    z = mu - stdev*log(yp)}
  else {
    z = mu - (stdev/sh)*(1-yp^-sh)}
}


## ----return level GEV-----------------------------------------------------------------------------------------------------------
y1max = zp(par_GEV$mle[1], par_GEV$mle[2], par_GEV$mle[3], 1) #mu, stdev, sh, year
y1max


## ----comparação com amostra, warning=FALSE, message=FALSE-----------------------------------------------------------------------
#agrupar por ano
table = ws %>% group_by(ANO) %>% summarise(n=n(), max=max(FF_KM)); table

#média de velocidade máxima
mean_max = mean(table$max); mean_max


## ----graphs GEV return level, message=FALSE-------------------------------------------------------------------------------------
anos = c(1,2,3,4,5,7,10,15,20,30,40,50,75,100,200,300)
niveis = c(0)
for (i in 1:length(anos)){
  niveis[i] = zp(par_GEV$mle[1], par_GEV$mle[2], par_GEV$mle[3], i)
}



sem_log = ggplot() + 
  geom_point(aes(anos,niveis, color=niveis)) + 
  geom_smooth(method='lm', color = 'red') + 
  ggtitle('Velocidade Vento Máxima por ano') + 
  labs(x='ano', y='velocidade max') +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_classic() + 
  theme(legend.position="top")

com_log = ggplot() + 
  geom_point(aes(log(anos),niveis, color=niveis)) + 
  geom_smooth(aes(log(anos),niveis), method='lm', color = 'red') + 
  ggtitle('Velocidade Vento Máxima p/ log(ano)') + 
  labs(x='log(ano)', y='velocidade max') +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_classic() + 
  theme(legend.position="top")

sem_log+com_log



## ----return level GPD, warning=FALSE, message=FALSE-----------------------------------------------------------------------------

#média de observações por ano
m_per_year = mean(table$n); m_per_year


xm = function(mu, stdev, sh, year, xdat, u){
  
  #número de observações
  m = year*m_per_year #computar transformação para n observações
  
  k = length(xdat[xdat > u])
  n = length(xdat)
  #proportion of excesses
  zeta_u = k/n
  
  if(round(sh,1)==0) {
    xm_ = u + stdev * log(m*zeta_u)}
  else {
    xm_ = mu + stdev/sh * (((m*zeta_u)^sh)-1)}
}


xm1max = xm(par_GPD$mu, par_GPD$mle[1], par_GPD$mle[2], 1, ws$FF_KM, 53)
xm1max



## ----graph GPD ML, warning=FALSE, message=FALSE---------------------------------------------------------------------------------
anos_gpd = c(1,2,3,4,5,7,10,15,20,30,40,50,75,100,200,300)
niveis_gpd = c(0)
for (i in 1:length(anos_gpd)){
  niveis_gpd[i] = xm(par_GPD$mu, par_GPD$mle[1], par_GPD$mle[2], i, ws$FF_KM, 53)
}



sem_log_gpd = ggplot() + 
  geom_point(aes(anos_gpd,niveis_gpd, color=niveis_gpd)) + 
  geom_smooth(method='lm', color = 'red') + 
  ggtitle('Velocidade Vento Máxima p/ ano') + 
  labs(x='ano', y='velocidade max') +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_classic() + 
  theme(legend.position="top")

com_log_gpd = ggplot() + 
  geom_point(aes(log(anos_gpd),niveis_gpd, color=niveis_gpd)) + 
  geom_smooth(aes(log(anos_gpd),niveis_gpd), method='lm', color = 'red') + 
  ggtitle('Velocidade Vento Máxima p/ log(ano)') + 
  labs(x='log(ano)', y='velocidade max') +
  scale_color_gradient(low = "yellow", high = "red") +
  theme_classic() + 
  theme(legend.position="top")

sem_log_gpd+com_log_gpd

