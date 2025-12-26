set.seed(76)

#-----Ejercicio 1-----
#1.5
# cantidad de muestras
ns <- c(5, 10, 15, 20, 30, 50, 100, 200)



Nrep <- 1000 # cantidad de veces que calculamos el intervalo para cada n
alfa <- 0.05
tita <- 0.25
zalfa2 <- qnorm(0.975)
hits_por_n <- c(0,0,0,0,0,0,0,0)
i <- 1

for(n in ns){ # iteramos por los n
  hits <- 0
  for(j in 1:Nrep){ # iteramos por las Nrep
    muestra <- rbinom(n,1,tita)
    xn <- mean(muestra)
    se <- sqrt(xn*(1-xn))
    lower <- xn - (zalfa2 * se /sqrt(n))
    upper <- xn + (zalfa2 * se /sqrt(n))
    if(tita >= lower & tita <= upper){
      hits <- hits + 1
    }
  }
  hits_por_n[i] <- hits
  i <- i + 1
}

cobertura_empirica <- hits_por_n / Nrep

#graficamos
plot(ns, cobertura_empirica,
     type = "b",
     pch = 19,
     col = "blue",
     xlab = "tamaño de la muestra",
     ylab = "cobertura empírica",
     xaxt = "n",
     main = "Cobertura empírica de los intervalos de confianza asintóticos
     según el tamaño de la muestra")

axis(1, at = ns, labels = ns)

#------ Ejercicio 2-----
set.seed(76)
# 2.3
Se <- 0.9
Sp <- 0.95
tita <- 0.25

# formula de p
# p <- Se * tita + (1 - Sp) * (1 - tita)

#a) dejando fijos Sp, Se, variamos tita
tita_vals <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # 9 vals

p_vals <- Se * tita_vals + (1-Sp) * (1-tita_vals)

plot(tita_vals, p_vals, axes=FALSE,
     type = "b",
     pch = 19,
     col = "blue",
     xlab = "tita",
     ylab = "p",
     main = "Valor de p en función de tita con Se = 0.9 y Sp = 0.95")
axis(1, at = tita_vals)


axis(2, at = p_vals)
box()

# b) variando Se

se_vals <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # 9 vals

p_vals <- se_vals * tita + (1-Sp) * (1-tita)

plot(se_vals, p_vals, axes=FALSE,
     type = "b",
     pch = 19,
     col = "blue",
     xlab = "Se",
     ylab = "p",
     main = "Valor de p en función de Se con tita = 0.25 y Sp = 0.95")
axis(1, at = se_vals)
axis(2, at = p_vals)
box()

# c) variando Sp

sp_vals <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9) # 9 vals



p_vals <- Se * tita + (1-sp_vals) * (1-tita)

plot(sp_vals, p_vals, axes=FALSE,
     type = "b",
     pch = 19,
     col = "blue",
     xlab = "Sp",
     ylab = "p",
     main = "Valor de p en función de Sp con tita = 0.25 y Se = 0.9")
axis(1, at = sp_vals)
axis(2, at = p_vals)
box()
#-----2.1-----
set.seed(76)
# 2.1.6
#consideramos los valores de Se, Sp y Theta fijados en el 
#ejercicio 2.3
Se <- 0.9
Sp <- 0.95
tita <- 0.25

#calculamos el ECM del estimador de momentos del test imperfecto 
#para distintos valores de n (el tamaño de la muestra)
p<- Se*tita +(1-Sp)*(1-tita)
ECMs_imp<-c(0,0,0,0,0,0,0,0)
i<-1
for (n in ns){ #calculamos el ECM por cada n en ns
  ECM<- (1/n)*((1/Se+Sp+1)**2)*(p*(1-p))
  ECMs_imp[i]<- ECM
  i<-i+1
  
  
}

#graficamos el ECM del estimador de momentos del test imperfecto 
#en función de n (el tamaño de la muestra)

plot(ns, ECMs_imp,
     type = "b",
     pch = 19,
     col = "blue",
     xlab = "tamaño de la muestra(n)",
     ylab = "ECM",
     main = "ECM del estimador de momentos para el test imperfecto 
     según el tamaño de la muestra")

#Comparamos el ECM del estimador de momentos del test imperfecto con el ECM del
#estimador de momentos del test perfecto 

ECMs_perf <- c(0,0,0,0,0,0,0,0)
i<-1
for (n in ns){ #calculamos el ECM por cada n en ns
  ECM<- (1/n)*tita*(1-tita)
  ECMs_perf[i]<- ECM
  i<-i+1
  
  
}
#graficamos ECMs_perf y ECMs_imp en función a n para comparar el comportamiento 
#de las curvas

#n vs ECM_imp
plot(ns, ECMs_imp,
     type = "b",
     lwd=2,
     pch = 19,
     col = "blue",
     xlab = "tamaño de la muestra(n)",
     ylab = "ECM",
     main = "Comparación ECM: test perfecto vs imperfecto")

#agregamos la curva n vs ECM_perf
lines(ns, ECMs_perf,type="b",col="red", lwd=2, lty=1,pch=19)


legend("topright", legend=c("Test perfecto", "Test imperfecto"),
       col=c("red","blue"), lwd=2)


# 2.1.7
#simulaciones para estimar los valores de sesgo, varianza y ECM 
#del estimador de momentos 
set.seed(76)
# función para el estimador de momentos
theta_MoM <- function(x) {
  (mean(x)-1+Sp )/(Se-1+Sp)
}

# guardar resultados
resultados <- data.frame(n=integer(), sesgo=double(), var=double(), ecm=double())

for (n in ns) { #iteramos por los valores de ns
  tita_est <- numeric(Nrep)
  
  for (i in 1:Nrep) {
    # simulamos muestra Bernoulli con parametro p
    x <- rbinom(n, size=1, prob=p)
    #estimador de p segun la muestra i-esima
    tita_est[i] <- theta_MoM(x) 
  }
  #calculamos sesgo, varianza y ECM empiricos con 
  #la muestra i-esima
  sesgo <- mean(tita_est) - tita
  var <- var(tita_est)
  ecm <- var + sesgo^2
  
  resultados <- rbind(resultados, data.frame(n=n, sesgo=sesgo, var=var, ecm=ecm))
}

print(resultados)

#graficamos las aproximaciones empiricas del sesgo, la var 
#y el ecm en funcion al tamaño de muestra 

plot(resultados$n, resultados$sesgo, 
     type="b", 
     lwd="2", 
     col="blue", 
     xlab="Tamaño de la muestra(n)",
     ylab="Sesgo estimado",
     main="Sesgo empírico del estimador del momentos para el test 
     imperfecto en función al tamaño de la muestra") 
plot(resultados$n, resultados$var, 
     type="b", 
     lwd="2", 
     col="blue", 
     xlab="Tamaño de la muestra(n)",
     ylab="Varianza estimada",
     main="Varianza empírica del estimador del momentos para el test 
     imperfecto en función al tamaño de la muestra")
plot(resultados$n, resultados$ecm, 
     type="b", 
     lwd="2", 
     col="blue", 
     xlab="Tamaño de la muestra(n)",
     ylab="ECM estimado",
     main="ECM empírico del estimador del momentos para el test 
     imperfecto en función al tamaño de la muestra")

# Obs: los graficos son coherentes con los resultados teoricos

# 2.1.8
set.seed(76)
#distribucion del estimador de momentos para el test imperfecto
#usando Bootstrap paramétrico
n<- 10
B<-1000 #cant de muestras bootstrap
tita_bootstrap<-numeric(B) #vector de estimadores
p_true<- Se*tita +(1-Sp)*(1-tita)
#simulamos la "muestra original" con el
#p verdadero
muestra0<- rbinom(n, size=1, prob=p_true) 
tita0_hat<-theta_MoM(muestra0) #estimador de theta
p0_hat<- Se*tita0_hat+(1-Sp)*(1-tita0_hat) #calculamos el estimador plug in de p la muestra0

for (i in 1:B){ #simulacion bootstrap
  muestra<-rbinom(n,size=1, prob=p0_hat) #muestra bootstrap
  tita_bootstrap[i] <- theta_MoM(muestra) #estimador bootstrap
}

#histograma de los valores estimados de p con las muestras bootstrap
hist(tita_bootstrap, col = "steelblue", border = "white",freq=FALSE,
     main = expression("Distribución bootstrap de " ~ hat(theta)[MoM]),
     xlab = expression(hat(theta)[MoM]),
     ylab="Densidad")
abline(v = tita0_hat, col = "red", lwd = 2)  # valor del estimador original

#El histograma presenta huecos debido a que, al tratarse de una muestra 
#Bernoulli de tamaño pequeño, n=10,
#el estimador solo puede tomar un conjunto finito de valores discretos

#veamos con n=1000 (simplemente para ver como evoluciona el histograma
#y si efectivamente va tendiendo a algo continuo)

n<- 1000
B<-1000 #cant de muestras bootstrap
tita_bootstrap<-numeric(B) #vector de estimadores
p_true<- Se*tita +(1-Sp)*(1-tita)
muestra0<- rbinom(n, size=1, prob=p_true) #simulamos la "muestra original" con el
#p verdadero
tita0_hat<-theta_MoM(muestra0) #estimador de theta
p0_hat<- Se*tita0_hat+(1-Sp)*(1-tita0_hat) #calculamos el estimador de momentos de p la muestra0

for (i in 1:B){ #simulacion bootstrap
  muestra<-rbinom(n,size=1, prob=p0_hat) #muestra bootstrap
  tita_bootstrap[i] <- theta_MoM(muestra) #estimador bootstrap
}

#histograma de los valores estimados de p con las muestras bootstrap (n=1000)
hist(tita_bootstrap, col = "steelblue", border = "white",freq=FALSE,
     main = expression("Distribución bootstrap de " ~ hat(theta)[MoM]),
     xlab = expression(hat(theta)[MoM]),
     ylab="Densidad")
abline(v = tita0_hat, col = "red", lwd = 2)  # valor del estimador original

#efectivamente, vemos que va tendiendo a lo que pareciera una normal. 

#-----Ejercicio2.2-----

# 2.2.9
set.seed(76)

tita <- 0.25
Se <- 0.9
Sp <- 0.95
p <- Se*tita+(1-Sp)*(1-tita)
n <- c(10,20,30,50,100,200) 

B <- Nrep
nfilas <- 6
intervalos_percentil <- data.frame(
  n = numeric(nfilas),
  lim_inf = numeric(nfilas),
  lim_sup = numeric(nfilas)
)
fila <- 1


for (i in n){
  #simulamos la "muestra original" con el p verdadero
  muestra_0 <- rbinom(i,1,p)
  tita_0 <- theta_MoM(muestra_0)
  p0<- tita_0 * Se + (1-Sp)*(1-tita_0)
  theta_MOM_bootstrap <- numeric(B)
  for (j in 1:B){
    muestra <- rbinom(i,1,p0)
    tita_actual <- theta_MoM(muestra)
    theta_MOM_bootstrap[j] <- tita_actual
  }
  lim.inf <- quantile(theta_MOM_bootstrap, 0.025)
  lim.sup <- quantile(theta_MOM_bootstrap, 0.975)
  intervalos_percentil$n[fila] <- i
  intervalos_percentil$lim_inf[fila] <- lim.inf
  intervalos_percentil$lim_sup[fila] <- lim.sup
  fila <- fila +1
} 

#graficamos los intervalos

library(ggplot2)

ggplot(intervalos_percentil, aes(y = n)) +
  geom_segment(aes(x = lim_inf, xend = lim_sup, y = n, yend = n), 
               color = "steelblue", size = 1) +
  geom_point(aes(x = lim_inf), color = "steelblue", size = 2) +
  geom_point(aes(x = lim_sup), color = "steelblue", size = 2) +
  geom_vline(xintercept = tita, color = "green", linetype = "dashed", size = 1) +
  labs(title = "Intervalos de Confianza Bootstrap por tamaño de muestra",
       subtitle = paste("Línea verde discontinua indica tita =", tita),
       x = "Límites del Intervalo",
       y = "Tamaño de muestra (n)") +
  theme_minimal()


# 2.2.10
#intervalos de confianza asintoticos 0.95
set.seed(76)
tita <- 0.25
Se <- 0.9
Sp <- 0.95
p <- Se*tita+(1-Sp)*(1-tita)
ns <- c(10,20,30,50,100,200) 


nfilas <- 6
intervalos_asintoticos <- data.frame(
  n = numeric(nfilas),
  lim_inf = numeric(nfilas),
  lim_sup = numeric(nfilas)
)
fila <- 1

for(n in ns){
  muestra <- rbinom(n,1,p)
  tita_est <- theta_MoM(muestra)
  t_line <- mean(muestra)
  zalfa2 <- qnorm(0.975)
  media_long <- (zalfa2 /(Se-1+Sp)) * sqrt(t_line * (1-t_line)) / sqrt(n)
  lim_inf <- tita_est - media_long
  lim_sup <- tita_est + media_long
  intervalos_asintoticos$n[fila] <- n
  intervalos_asintoticos$lim_inf[fila] <- lim_inf
  intervalos_asintoticos$lim_sup[fila] <- lim_sup
  fila <- fila +1
}


#graficamos
ggplot(intervalos_asintoticos, aes(y = n)) +
  geom_segment(aes(x = lim_inf, xend = lim_sup, y = n, yend = n), 
               color = "steelblue", size = 1) +
  geom_point(aes(x = lim_inf), color = "steelblue", size = 2) +
  geom_point(aes(x = lim_sup), color = "steelblue", size = 2) +
  geom_vline(xintercept = tita, color = "green", linetype = "dashed", size = 1) +
  labs(title = "Intervalos de Confianza asintóticos por tamaño de muestra",
       subtitle = paste("Línea verde discontinua indica tita =", tita),
       x = "Límites del Intervalo",
       y = "Tamaño de muestra (n)") +
  theme_minimal()



# 2.2.11
#cuantos intervalos de cada tipo vamos a calcular para cada n
simulaciones <- 1000
set.seed(76)
tita <- 0.25
Se <- 0.9
Sp <- 0.95
p <- Se*tita+(1-Sp)*(1-tita)
ns <- c(10,20,30,50,100,200) 

#primero para asintoticos
nfilas <- 6

intervalos_asintoticos <- data.frame(
  n = numeric(nfilas),
  longitud_promedio = numeric(nfilas),
  hits_promedio = numeric(nfilas)
)
fila <- 1

for(n in ns){
  hits <- 0
  longitudes <- 0
  for(i in 1:simulaciones){
    muestra <- rbinom(n,1,p)
    tita_est <- theta_MoM(muestra)
    t_line <- mean(muestra)
    zalfa2 <- qnorm(0.975)
    media_long <- zalfa2 /(Se-1+Sp) * sqrt(t_line * (1-t_line)) / sqrt(n)
    lim_inf <- tita_est - media_long
    lim_sup <- tita_est + media_long
    hit <- if(tita>=lim_inf & tita<= lim_sup){1}else{0}
    longitud <- lim_sup-lim_inf
    hits<- hits + hit
    longitudes <- longitudes + longitud
  }
  hits <- hits/simulaciones
  longitudes <- longitudes/simulaciones
  intervalos_asintoticos$n[fila] <- n
  intervalos_asintoticos$longitud_promedio[fila] <- longitudes
  intervalos_asintoticos$hits_promedio[fila] <- hits
  fila <- fila +1
}

#ahora para boostrap
intervalos_percentil <- data.frame(
  n = numeric(nfilas),
  longitud_promedio = numeric(nfilas),
  hits_promedio = numeric(nfilas)
)
fila <- 1


B <- Nrep


for(n in ns){
  hits <- 0
  longitudes <- 0
  for(i in 1:simulaciones){
    muestra0<- rbinom(n, size=1, prob=p) #simulamos la "muestra original" con el p verdadero
    tita_0 <- theta_MoM(muestra0)
    p0<- tita_0 * Se + (1-Sp)*(1-tita_0)
    theta_MOM_bootstrap <- numeric(B)
    for (j in 1:B){
      muestra <- rbinom(n,1,p0)
      tita_actual <- theta_MoM(muestra)
      theta_MOM_bootstrap[j] <- tita_actual
    }
    lim_inf <- quantile(theta_MOM_bootstrap, 0.025)
    lim_sup <- quantile(theta_MOM_bootstrap, 0.975)
    hit <- if(tita>=lim_inf & tita<= lim_sup){1}else{0}
    longitud <- lim_sup-lim_inf
    hits<- hits + hit
    longitudes <- longitudes + longitud
  }
  hits <- hits/simulaciones
  longitudes <- longitudes/simulaciones
  intervalos_percentil$n[fila] <- n
  intervalos_percentil$longitud_promedio[fila] <- longitudes
  intervalos_percentil$hits_promedio[fila] <- hits
  fila <- fila +1
}

intervalos_asintoticos$metodo <- "Asintotico"
intervalos_percentil$metodo <- "Percentil"
datos_combinados <- rbind(intervalos_asintoticos, intervalos_percentil)
#graficamos
ggplot(datos_combinados, aes(x = n, y = longitud_promedio, color = metodo)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Comparación de Longitud Promedio",
    subtitle = "Intervalos Asintóticos vs. Percentil",
    x = "Tamaño de muestra (n)",
    y = "Longitud Promedio",
    color = "Método"
  ) +
  theme_minimal()


ggplot(datos_combinados, aes(x = n, y = hits_promedio, color = metodo)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  ylim(0, 1) + 
  labs(
    title = "Comparación de Hits Promedio (Cobertura)",
    subtitle = "Intervalos Asintóticos vs. Percentil",
    x = "Tamaño de muestra (n)",
    y = "Hits Promedio",
    color = "Método"
  ) +
  theme_minimal()

#-----Ejercicio 2.3-----

# 2.3.13

#simulaciones para estimar los valores de sesgo, varianza y ECM del estimador 
#de momentos truncado 

# función para estimador de momentos truncado
theta_MoM <- function(x) {
  (mean(x)-1+Sp )/(Se-1+Sp)
}

theta_MoM_truncado <- function(x){
  if (theta_MoM(x) < 0){
    return(0)
  }
  else if (theta_MoM(x)> 1){
    return(1)
  } 
  else { 
    return(theta_MoM(x))
  }
}

# Aproximamos con simulaciones el sesgo, la varianza y el ECM para distintos valores de n
n <- c(10,20,30,100,400,800,1000)
resultados <- data.frame(n=integer(), sesgo=double(), var=double(), ecm=double())

set.seed(76) # nos quedamos con el primer resultado de esta semilla
for (i in n) { #iteramos por los valores de n
  p_est <- numeric(Nrep)
  
  for (j in 1:Nrep) {
    # simulamos muestra Bernoulli con parametro p
    x <- rbinom(n, size=1, prob=p)
    p_est[j] <- theta_MoM_truncado(x) #estimador de theta segun la muestra i-esima
  }
  #calculamos sesgo, varianza y ECM con la muestra i-esima
  sesgo <- mean(p_est) - p
  var <- var(p_est)
  ecm <- var + sesgo^2
  
  resultados <- rbind(resultados, data.frame(n=i, sesgo=sesgo, var=var, ecm=ecm))
}

print(resultados)


#graficamos las aproximaciones empiricas del sesgo, la var y el ecm en 
#funcion al tamaño de muestra 

plot(resultados$n, resultados$sesgo, 
     type="b", 
     lwd="2", 
     col="blue", 
     xlab="n",
     ylab="sesgo estimado",
     ylim = c(-0.08, 0.08),
     main="Sesgo del estimador de momentos truncado para distintos tamaños
     de la muestra")

plot(resultados$n, resultados$var, 
     type="b", 
     lwd="2", 
     col="blue", 
     xlab="n",
     ylab="var estimada",
     main="Varianza del estimador de momentos truncado para distintos tamaños
     de la muestra")

plot(resultados$n, resultados$ecm, 
     type="b", 
     lwd="2", 
     col="blue", 
     xlab="n",
     ylab="ECM estimado",
     main="ECM del estimador de momentos truncado para distintos tamaños
     de la muestra")


#obs no sabemos exacto si es insesgado , vemos si es asint insesgado con valores mas grandes de n 
n <- c(1000,2000,5000,10000)
resultados <- data.frame(n=integer(), sesgo=double(), var=double(), ecm=double())

set.seed(76) # nos quedamos con el primer resultado de esta semilla
for (i in n) { #iteramos por los valores de n
  p_est <- numeric(Nrep)
  
  for (j in 1:Nrep) {
    # simulamos muestra Bernoulli con parametro p
    x <- rbinom(n, size=1, prob=p)
    p_est[j] <- theta_MoM_truncado(x) #estimador de theta segun la muestra i-esima
  }
  #calculamos sesgo, varianza y ECM con la muestra i-esima
  sesgo <- mean(p_est) - p
  var <- var(p_est)
  ecm <- var + sesgo^2
  
  resultados <- rbind(resultados, data.frame(n=i, sesgo=sesgo, var=var, ecm=ecm))
}

print(resultados)
#es asint insesgado


#distribucion asintotica con estimacion a traves de simulaciones Monte Carlo

set.seed(76) 
n <- c(10,20,30,100,400,800,1000,5000,10000)
for (i in n){
  p_vector <-numeric(Nrep) #vector de estimadores
  for (j in 1:Nrep){ #simulacion 
    muestra<-rbinom(n,size=1, prob=p) #muestra aleatoria
    p_vector[j] <- theta_MoM_truncado(muestra) #estimador 
    
  }
  titulo_dinamico <- paste0("Distribucion de theta_MoM_truncado para n=", i)
  #histograma de los valores estimados de p con las muestras aleatorias
  hist(p_vector, breaks =21, col = "steelblue", border = "white", freq=TRUE,
       main = titulo_dinamico,
       xlab = expression(hat(Tetha_MoM_truncado)))
  plot(density(p_vector, bw=0.1),
       main = titulo_dinamico,
       xlab = expression(hat(Tetha_MoM_truncado)))
}

#------Ejercicio 3-----
set.seed(76)
Se<-0.9
Sp<-0.95
n_pre<-100
n_post<-100 
theta_pre<-0.2
theta_post<-0.15
p_pre<-(Se+Sp-1)*theta_pre+(1-Sp)
p_post<-(Se+Sp-1)*theta_post+(1-Sp)
alpha<-0.05
z_alpha2<-qnorm(0.975)

#pivote del test 
pivote <- function(x1, x2) {
  #tamaños muestrales
  npre  <- length(x1)
  npost <- length(x2)
  
  # estimadores de momentos
  tita_mom_pre  <- theta_MoM(x1)
  tita_mom_post <- theta_MoM(x2)
  
  dif <- tita_mom_post - tita_mom_pre
  
  # estinadores de p_pre y p_post
  p_pre_hat  <- mean(x1)
  p_post_hat <- mean(x2)
  
  # var estimada
  sigma_hat <- (p_pre_hat * (1 - p_pre_hat)) / (npre  * (Se + Sp - 1)^2) +
    (p_post_hat * (1 - p_post_hat)) / (npost * (Se + Sp - 1)^2)
  
  # estadístico pivote
  res <- dif / sqrt(sigma_hat)
  
  return(res)
}


#definimos el test de hipotesis
test <- function(x1,x2){
  #valor del estadítico observado
  t_obs <- pivote(x1,x2) 
  
  # no rechazo si no está bien definido
  if (is.na(t_obs) || is.infinite(t_obs)) {
    return(0)  
  }
  
  #condición de rechazo
  if (abs(t_obs) > z_alpha2){
    return(1)
    
  } else {
    return(0)
  }
}


#aplicamos el test a distintos tamaños de muestra
ns<-c(10,20,30,50,100,1000,5000,10000)
Nreps<- 1000
prop_rechazos<-c(0,0,0,0,0,0)
j<-1
for (n in ns){
  contador<-0
  for (i in (1:Nreps)){
    x1<- rbinom(n,1,p_pre )
    x2<- rbinom(n,1,p_post )
    contador<- contador+test(x1,x2)
    
  }
  prop_rechazos[j]<-contador/Nreps
  j<-j+1
}
print(prop_rechazos)

plot(ns, prop_rechazos,type = "b", col = "red",
     xlab="Tamaño de la muestra(n)", ylab="Proporción de rechazos",
     main="Proporción de rechazos del test asintótico bajo H1")

#Intervalos de confianza 
#construcción y gráfico de intervalos.

# vector para guardar intervalos
intervalos <- data.frame(replica = integer(),
                         lim_inf = numeric(),
                         lim_sup = numeric())

for (i in 1:100) {
  muestra1 <- rbinom(n_pre, 1, p_pre)
  muestra2 <- rbinom(n_post, 1, p_post)
  
  tita_mom_1 <- theta_MoM(muestra1)
  tita_mom_2 <- theta_MoM(muestra2)
  dif <- tita_mom_2 - tita_mom_1
  
  p_pre_hat  <- mean(muestra1)
  p_post_hat <- mean(muestra2)
  
  var_muestra1 <- (p_pre_hat * (1 - p_pre_hat)) / (n_pre  * (Se + Sp - 1)^2)
  var_muestra2 <- (p_post_hat * (1 - p_post_hat)) / (n_post * (Se + Sp - 1)^2)
  
  se_delta <- sqrt(var_muestra1 + var_muestra2)
  
  lim_inf <- dif - z_alpha2 * se_delta
  lim_sup <- dif + z_alpha2 * se_delta
  
  intervalos <- rbind(intervalos,
                      data.frame(replica = i,
                                 lim_inf = lim_inf,
                                 lim_sup = lim_sup))
}

# valor verdadero de Delta
Delta_true <- theta_post-theta_pre

# gráfico
library(ggplot2)
ggplot(intervalos, aes(x = replica)) +
  geom_errorbar(aes(ymin = lim_inf, ymax = lim_sup), width = 0.2, color = "steelblue") +
  geom_hline(yintercept = Delta_true, color = "red", linetype = "dashed", size = 1) +
  labs(title = "Intervalos de confianza simulados para Δ",
       x = "Simulación",
       y = "Intervalo de confianza") +
  theme_minimal()



# Cobertura empírica del test
#Para simular nivel empirico, consideramos un caso bajo H0
#Consideramos theta_pre=theta_post=0.2-> p_pre=p_post=0.62 
#Seguimos considerando Sp y Se como definimos antes

#Tamaños de muestra distintos 
ns_pre  <- c(5,10,20,30,50,60,70,100,500,1000,1500, 2000,3000,3500, 5000,6000,8000,10000)
ns_post <- c(3,7,15,35,47,61,65,98,503,995,1450, 1997,3000,3504, 4900,5010,8000,9999)
p_pre<-0.62
p_post<-0.62
prop_rechazos <- numeric(length(ns_pre))

for (k in seq_along(ns_pre)) {
  npre  <- ns_pre[k]
  npost <- ns_post[k]
  
  rechazos <- 0
  for (i in 1:Nreps) {
    muestra1 <- rbinom(npre, 1, p_pre)
    muestra2 <- rbinom(npost, 1, p_post)
    
    rechazos <- rechazos + test(muestra1, muestra2)
  }
  
  # proporción de rechazos
  prop_rechazos[k] <- rechazos / Nreps
}

prop_rechazos


#grafico 

df <- data.frame(
  npre = ns_pre,
  npost = ns_post,
  nivel = prop_rechazos
)

ggplot(df, aes(x = npre, y = npost, color = nivel, size = nivel)) +
  geom_point(alpha = 0.9) +
  scale_color_gradient2(
    low = "pink",   # violeta fuerte
    mid = "orange",     # blanco en 0.05
    high = "#B22222",  # rojo intenso
    midpoint = 0.05
  ) +
  scale_size(range = c(3, 8)) +
  labs(x = "n pre", y = "n post", color = "Nivel empírico", size = "Nivel empírico",
       title = "Nivel empírico en función de n") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey80"),  # grilla principal gris
    panel.grid.minor = element_line(color = "grey90"),  # grilla secundaria más clara
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


#hacemos lo mismo pero considerando las dos muetsras de tamaño igual

ns_pre  <- c(5,10,20,30,50,60,70,100,500,1000,1500, 2000,3000,3500, 5000,6000,8000,10000)
ns_post <- c(5,10,20,30,50,60,70,100,500,1000,1500, 2000,3000,3500, 5000,6000,8000,10000)
p_pre<-0.62
p_post<-0.62
prop_rechazos <- numeric(length(ns_pre))

for (k in seq_along(ns_pre)) {
  npre  <- ns_pre[k]
  npost <- ns_post[k]
  
  rechazos <- 0
  for (i in 1:Nreps) {
    muestra1 <- rbinom(npre, 1, p_pre)
    muestra2 <- rbinom(npost, 1, p_post)
    
    rechazos <- rechazos + test(muestra1, muestra2)
  }
  
  # proporción de rechazos
  prop_rechazos[k] <- rechazos / Nreps
}

prop_rechazos


#cuando los tamaños de las muestras son iguales estamos mas cerca del 0.05. Cuando
#los tamaños son distintos, el test es más "conservador".



#graficamos el nivel empirico en función a los tamaños de la muestra
plot(ns_pre, prop_rechazos, col="blue", main="Nivel empírico del test 
     asintótico para distintos valores de muestras",
     xlab="Tamaño de las muestras",ylab="Nivel empírico" )

