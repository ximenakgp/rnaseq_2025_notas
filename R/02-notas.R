# Bioconductor
## Fundador Robert Gentleman
## Todos los packages de Bioconductor tienen una vignette
## Los paquetes de bioconductor los evaluan de forma diaria y los de cran solo cuando se suben
## Release: es como la versión para consumo público
## Devel: Versión en desarrollo (pero es es más para desarrolladores)

BiocManager::version()
BiocManager::valid()

#Cuando seleccionamos un paquete y lo abrimos en help podemos ver su documentacion
?reconut3::create_rse # El signo también funciona cuando queremos acceder a la documentacion de la funcion en ese paquete

.libPaths()
# Nos dice en donde esta instalado el paquete
# [1] "/home/ximenagp/R/x86_64-pc-linux-gnu-library/4.4"
# [2] "/usr/local/lib/R/site-library"
# [3] "/usr/lib/R/site-library"
# [4] "/usr/lib/R/library"
# En algunas direcciones tenemos permiso para instalar y en otras no

?BiocManager::install
?install.packages
# El argumento lib no dice en que lugar se estan instalando

# Cada paquete es un directorio

BiocManager::install()
length(dir(.libPaths()))
#[1] 541 Paquetes actuales


# Cada 6 meses la version devel se vuelve la release y le cambian el numero a la version
