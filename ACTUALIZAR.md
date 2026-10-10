# Actualización para StatatoR en GitHub

1. Descomprimir el ZIP.
2. Subir su contenido a la raíz del repositorio, conservando la estructura:
   - DESCRIPTION (dependencias ivreg/fixest; versión 0.3.0).
   - NAMESPACE (los tres métodos S3 nuevos, importaciones y exportaciones previas).
   - R/ext_summary.R (reemplaza la función anterior completa).
   - R/package_imports.R (metadatos de importaciones para roxygen).
   - man/ (ayuda generada, incluido ext_summary.Rd).
3. Guardar el commit en GitHub y reiniciar R antes de actualizar la biblioteca.

```r
remotes::install_github('JISXII/StatatoR',dependencies=TRUE,force=TRUE)
library(statator)
help('ext_summary',package='statator')
```

Si no se usa el paquete, basta source de R/ext_summary.R; las dependencias
que requiera cada modelo deben estar instaladas.

La actualización preserva las funciones existentes box.cox, ladder, gladder
y qladder. No modifica el repositorio remoto automáticamente.

La declaración MIT del DESCRIPTION original sigue igual. R CMD check
completa sin errores ni warnings; mantiene una nota sobre añadir LICENSE.

La versión 0.3.0 incorpora margins(), marginsplot(), iia_test() y nlcom().
Subir TODOS los archivos R y de ayuda del ZIP: las nuevas funciones comparten utilidades privadas.
README.md explica el alcance, los argumentos y ejemplos de uso.
