# Márgenes, gráficos y tests en StatatoR 0.3.0

Estas funciones son independientes de los nombres y las categorías de la clase.
Su generalidad está en los inputs y en el cálculo a partir del modelo, no en
suponer que cualquier procedimiento estadístico sirve para cualquier estimador.
Cada función documenta su alcance y rechaza las combinaciones que no admite.

## Instalar la versión local

El ZIP actualizado todavía no se ha publicado en GitHub. Para probarlo ahora:

```r
# install.packages("remotes")  # solo si falta
remotes::install_local(file.choose(), dependencies = TRUE, force = TRUE)
```

Seleccionar `StatatoR-update.zip`. Reiniciar R si había una versión anterior
cargada y ejecutar:

```r
library(statator)
help(margins, package = "statator")
```

Después de subir el contenido del ZIP a la raíz del repositorio, puede usarse
la instalación habitual desde `JISXII/StatatoR`. No basta con subir un solo R:
las funciones comparten utilidades privadas y requieren los nuevos NAMESPACE,
DESCRIPTION y archivos de ayuda.

## Uso habitual: una o dos líneas por análisis

```r
margins(modelo)                                  # predicción promedio
margins(modelo, atmeans = TRUE)                  # predicción en las medias
margins(modelo, dydx = "edad")                   # efecto promedio
margins(modelo, dydx = "edad", atmeans = TRUE)   # efecto en las medias
margins(modelo, change = list(tratamiento = c(0, 1)))
margins(modelo, at = list(edad = c(20, 40, 60)))
```

`at` fija los valores indicados y conserva la distribución observada de las
otras variables. Con `atmeans = TRUE`, las otras numéricas se fijan en sus medias
y las categóricas en sus modas. Las binarias numéricas también se fijan en su
media; no se interpreta su media como un individuo realmente observado.

`average = FALSE` devuelve resultados por observación/escenario. `by` permite
promediar por grupos. `newdata` recibe datos explícitos para diseños más complejos.
`outcome` selecciona las categorías de respuesta. Varias variables en `at`
generan todas las combinaciones, como una grilla de escenarios.

Para una variable categórica, declarar los dos valores del cambio es lo más
transparente. En modelos ordinarios, `dydx` sigue las convenciones de
marginaleffects, que trata las binarias/categóricas mediante cambios discretos.
En el motor de elección, `dydx` deriva las numéricas, incluidas las binarias
numéricas, y `change` pide explícitamente el cambio discreto.

## Covarianza, precisión e inferencia

```r
resumen <- ext_summary(modelo, vce = "robust")
margins(resumen, dydx = "edad")
```

El objeto conserva el modelo, su covarianza y la distribución de referencia.
También puede indicarse directamente `vce`, `hc_type`, `cluster`, `correction`
o `vcov`, como en ext_summary. `vcov` y las opciones VCE son alternativas.
Para modelos fuera del alcance de ext_summary se usa la covarianza nativa
o una solicitud explícita admitida por marginaleffects.

Los promedios ordinarios dan el mismo peso a cada fila por defecto. Si el modelo
usa ponderadores y se busca un promedio ponderado, indicar `wts = TRUE` o
entregar el vector apropiado mediante `wts`; son argumentos de marginaleffects
que se reenvían explícitamente. No se asume que cualquier promedio ponderado
deba tener la misma interpretación que un margen de encuesta en Stata.

`level` admite 0.90 o 90. `df = Inf` usa z y la normal; un df finito usa t.
`digits` afecta la impresión, no los cálculos. `numderiv = "richardson"` es el
predeterminado; `eps` controla el paso de diferenciación de las variables.
Una escala muy distinta puede requerir comprobar sensibilidad a otros pasos.

En el ordinal de la clase se verificó:

```r
margins(m_ologit, dydx = "gpa", atmeans = TRUE, eps = 0.001)
```

Se usa la matriz completa de pendientes y cortes. En NB2 se conserva la
covarianza conjunta para nlcom, mientras que las predicciones utilizan el
bloque de las pendientes: la probabilidad/media predicha no depende directamente
del parámetro auxiliar de dispersión.

## Presentación y extracción

La estimación se identifica como `Probability`, `Prediction`, `dy/dx`, `Change`
o la elasticidad pedida (`ey/ex`, `ey/dx`, `dy/ex`). Las otras columnas son
`Std. Error`, `z`/`t`, `P>|z|`/`P>|t|`, `Sig.`, `CI Lower` y `CI Upper`.
Las categorías figuran como `Outcome`; los bloques se separan con espacios.

```r
r <- margins(modelo, at = list(edad = c(20, 40, 60)))
r$table              # tabla sin redondear
as.data.frame(r)     # la misma tabla
marginsplot(r, xvar = "edad")
```

El gráfico consume las estimaciones y los intervalos guardados, sin reestimar.
`group`, `facet`, `ci`, `connect`, `horizontal` y `title` permiten organizarlo.
Devuelve un ggplot que admite capas y temas adicionales.

## Modelos ordinarios y de elección

| Diseño | Cálculo y alcance |
|---|---|
| lm, GLM logit/probit, Poisson/NB2, multinom, polr, ivreg y fixest | Interfaz de marginaleffects; comprobados ejemplos de las clases principales. Otros modelos dependen del soporte del backend. |
| CLOGIT `type="pc1"` | Predicción normalizada dentro del estrato, con una elección por estrato. Solo predicciones; no medias que destruyan el conjunto de elección. |
| CLOGIT `type="pu0"` | Inversa logística del índice con efecto fijo cero; predicciones, derivadas y cambios discretos. Es una cantidad diferente de pc1. |
| MLOGIT ordinario usado como ASCLOGIT | Probabilidades, cambios discretos y efectos propios/cruzados; variables e interacciones obtenidas del modelo. |
| MLOGIT anidado | Predicciones y diferencias/derivadas numéricas, condicionadas a un ajuste identificado y covarianza válida. La comprobación numérica de la clase cubre su LR, no certifica los márgenes de su ajuste inestable. |

El motor CLOGIT/MLOGIT requiere modelos sin pesos ni offsets; MLOGIT debe
ser no restringido, determinista y tener índices de caso y alternativa.
No admite modelos mixtos, probit de elección, heteroscedásticos, panel o ranking.
Los datos de predicción MLOGIT deben incluir todas las alternativas en cada
escenario. Si el ajuste usa conjuntos desbalanceados, preparar un `newdata`
balanceado explícito. Los márgenes con efectos fijos requieren escenarios
válidos para los niveles de los efectos absorbidos; una media numérica de
un identificador no es un grupo válido.

Para ASCLOGIT se infiere si una variable cambia dentro del caso o es común
a la persona. Se puede declarar explícitamente:

```r
margins(modelo, dydx = "precio", atmeans = TRUE,
        variable_type = c(precio = "alternative"))
```

`Alternative` identifica el atributo cambiado y `Outcome` la probabilidad
afectada. Las medias de atributos se calculan por alternativa; las características
del caso se promedian una vez por caso. Para escenarios particulares, factores
complejos o una interpretación distinta de las medias, usar `newdata`.

## IIA

```r
iia_test(m_base)
iia_test(m_base, method = "hausman")
iia_test(m_base, omit = c("Prepaid", "Uninsure"))
```

Para multinom: Hausman-McFadden, Hausman con covarianza conjunta SUEST y
Small-Hsiao. Admite cualquier cantidad de categorías nominales desde tres,
factores, interacciones y transformaciones de regresores, manteniendo la
matriz de diseño estimada. Requiere `model=TRUE`, modelo sin ponderación,
penalización, offsets o máscaras personalizadas y diseños de rango completo.
Si una mitad pierde categorías, aparece separación/no convergencia o hay
singularidad, debe revisarse el diseño o la partición; no se fabrica un p-valor.

```r
i <- iia_test(m_base, seed = 12345)
iia_test(m_base, split = i$split)
i$table
```

El número de semilla no reproduce por sí solo la partición de Stata. El
estadístico Small-Hsiao de la clase difiere por ese motivo, como ya se documentó.
Los estadísticos de los tres tests coinciden con la implementación revisada de R.
La función marca diferencias de covarianzas Hausman negativas **o indefinidas**
como no interpretables, y omite su p-valor. En la clase, también la pequeña
Chi2 positiva al omitir Prepaid usa una matriz indefinida: se conserva Chi2,
pero no se muestra el p=1 nominal del PDF.

Para MLOGIT ordinario, se admite Hausman con omisión automática:

```r
iia_test(m_asc, terms = c("p", "q"), order = "full-reduced")
```

`terms` y el orden anterior reproducen la comparación del PDF; el orden
convencional predeterminado es `reduced-full`. SUEST/Small-Hsiao para MLOGIT
no están implementados y se rechazan expresamente. Para CLOGIT, entregar
el ajuste restringido ya estimado, con la misma especificación:

```r
iia_test(m_clogit, reduced_model = m_choc_rest, order = "full-reduced")
iia_test(m_anidado, method = "lr", reduced_model = m_asc)
```

El LR comprueba misma muestra/especificación y contrasta disimilitudes iguales
a uno. No corrige un ajuste anidado inadmisible como el del PDF. Una prueba
no significativa no demuestra que IIA sea verdadera.

## Combinaciones no lineales generales

```r
nlcom(modelo, function(b) c(OR_edad = exp(b["edad"])), null = 1)
```

La expresión puede devolver un vector de estimaciones diferenciables. También
se admiten un vector de parámetros con nombres y su matriz de covarianzas.
La salida incluye `table`, `vcov` de las estimaciones y `jacobian`.
Los intervalos Wald se calculan directamente en la escala transformada y
pueden salir de sus límites naturales. La prueba de nlcom sobre exp(b)=1
linealiza esa transformación; no es idéntica a la prueba sobre b=0 que
conserva ext_summary con eform=TRUE.

## Verificación y fuentes

Se aprobaron 60 comprobaciones iniciales: fórmulas analíticas, propagación de
covarianzas, inferencias de la clase, otros nombres y números de alternativas,
factores, interacciones, transformaciones, extracción y datos del gráfico.
El detalle está en `VALIDATION_POSTESTIMATION.csv` y en el script de validación.

- [marginaleffects: predicciones](https://marginaleffects.com/man/r/predictions.html).
- [marginaleffects: comparaciones](https://marginaleffects.com/man/r/comparisons.html).
- [marginaleffects: derivadas](https://marginaleffects.com/man/r/slopes.html).
- [Stata: CLOGIT postestimation](https://www.stata.com/manuals/rclogitpostestimation.pdf).
- [Long y Freese: tests IIA, STB-58](https://www.stata.com/products/stb/journals/stb58.pdf).
- [Stata: SUEST](https://www.stata.com/manuals/rsuest.pdf).
