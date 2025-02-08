## Creamos este archivo con el siguiente comando:
# usethis::use_r("01-notas.R")

usethis::create_github_token()
usethis::edit_r_environ() ## Linux
gitcreds::gitcreds_set() ## macOS o winOS

## Configura tu usuario de GitHub
usethis::edit_git_config()


## Queremos usar git con nuestro Rproj
usethis::use_git()

## Para conectar tu repositorio local de Git con los servidores de GitHub
usethis::use_github()
