# Metagenomics  
This project was developed by for academic use in the Core Lab Genomics investigation group at the Instituto Tecnológico y de Estudios Superiores de Monterrey (ITESM).

This program was designed for making a metagenomic analysis from 16s and shotgun sequencing. It is made from TXT or CSV format data from a taxonomic analysis from [KRAKEN 2](https://github.com/DerrickWood/kraken2.).

This program aims to create statiscal analyses and allow the user the visualization of taxonomical data from KRAKEN 2. 

# Libraries and Utils 

The firts part of the code is designed to declare different functions that will be used for the creation of barplots and boxplots that will be part of the results and analysis. 

It's important to run these first lines because some functions are not part of any library and they are necessary for the program to work correctly.

Aside from that, some libraries are installed too, these will be used to make different statistical anlysis. 


On the line:
```Rscript
geom_boxplot(fill=sample(mypal, length(treatcol)),fatten=1, outlier.shape = NA)
```
Random colors will be generated for the boxplots from Shannon and Simpson indexes. It is possible to change the parameter `fill =` in order to assign specific colors depending on the number of variables. This must be done following this format `fill = c(C1,C2,C3,...`. If the user wants only one color for the boxplots, it must be done like this  `fill = *hexadexcimal color code*`. If the user wants to have blank boxplots just needs to eliminate the `fill =` portion of the code. This will leave the code looking like this.

```Rscript
geom_boxplot(fill=sample(mypal, length(treatcol)),fatten=1, outlier.shape = NA)
```

It is recommended to use the hexadecimal code for the colors in order to get more specific and personalized colors depending on the study. If specificity is not required, the names of the colors can be used instead. 
### Nota 
Es muy importante mencionar que al modificar estos parámetros se tome en cuenta 

# Datos
En la segunda sección del código se declara el pathway de la carpeta donde se encuentran los datos con los cuales va a trabajar el programa. Es importante mencionar que dentro de esta carpeta se van a guardar los diferentes documentos e imágenes resultantes de este análisis. Por esto se recomienda destinar Carpetas específicas para cada análisis indívidual. 

Dentro de esta sección se hace un tratado de los datos para asegurar que el formato es el correcto y R o Rstudio no arroje códigos de error y códigos de advertencia. El input de este código debe de ser un archivo **.csv** o **.txt** proveniente del análisis taxonómico de [KRAKEN 2](https://github.com/DerrickWood/kraken2.). 

El input se debe de ver acomodado de la siguiente manera: 


| Rank | TaxId | Scientific Name | Sample 1 | Sample 2 | Sample 3 | 
| --- | --- | --- | --- | --- | --- |
| unclassified  | 0 | Uknown | 68 | 1232 | 1696 |
| superkingdom | 2 | Bacteria (bacteria) | 66708 | 76666 | 64937 | 
| genus | 131079 | Limnobacter| 0 | 0 | 0 |
| species | 2060312 | Altererythrobacter sp. B11 | 0 | 0 | 2 |

Es recomendable modificar solamente los nombres de las muestras (samples) para evitar cualquier inconveniente con el programa. Los nombres de las muestras pueden ser modificados libremente dependiendo de lo que se busque observar en el estudio. Más información de como se tratan, filtran y agrupan las muestras puede ser observado en la sección de [**Filtrado**](#Filtrado)
  
En la línea:
```Rscript
raw_data <- read.table(file.choose(), header = T, sep = "\t",quote = "", stringsAsFactors = F, fill = F)
```

Se debe de modificar el parámetro `"sep ="` dependiendo del tipo de documento que se esté usando: para **.txt** se usa `"\t"` y para **.csv** se usa `","`.

Las siguientes líneas de este código están dirigidas para asegurar que los datos se hayan introducido de manera correcta y no haya una pérdida de datos o se conviertan valores numéricos a strings y evitar que existan NA dentro de la tabla con la que se trabajará. 

La línea: 
```Rscript
raw_data <- raw_data[,c(n1,n2,n3)]
```
sirve para eliminar columnas que se hayan agregado en los pasos anteriores, ya sea al momento de cargar los datos, o al convertir los datos a números. Esta línea se puede modificar dependiendo del número de columnas que se hayan agregado. Ejemplo: si se agregaron dos columnas al final de un data frame que contenía 15 columnas originalmente, se tendrían que eliminar las columnas 16 y 17 (las dos que se agregaron). La línea quedaría de la siguiente manera `raw_data <- raw_data[,-c(16,17)]` de esta forma se eliminan las dos columnas agregadas por coerción. Esto no significa que haya una pérdida de datos, a veces las tablas que son arrojadas por KRAKEN contienen celdas invisibles que R las toma como columnas vacias y NA. 

Las líneas:
  
```Rscript
raw_data
str(raw_data)
```

funcionan para poder visualizar como fueron cargados los datos a R y ver si es que existe algún problema con los mismos. Si todo fue cargado de manera correcta se deberían de ver los datos de la misma manera que se ve el documento de origen. En el caso de `srt(raw_data)` se deben de ver los datos de la siguiente manera: 
            
```
 $ Rank           : chr  "unclassified" "superkingdom" "superkingdom" "superkingdom" ...
 $ TaxId          : int  0 2 10239 2157 2759 200783 200795 200918 28890 200930 ...
 $ Scientific.Name: chr  "Unknown" "Bacteria <bacteria>" "Viruses" "Archaea" ...
 $ Sample 1       : int  33369 25062 0 1 0 0 1 1 1 0 ...
 $ Sample 2       : int  184966 73818 0 1 0 1 39 20 1 3 ...
 $ Sample 3       : int  329303 75747 0 0 0 0 29 0 0 0 ...
```
Las siguientes líneas están destinadas a la creación de listas vacías para poder hacer el siguiente paso que es el filtrado de los datos a partir de 4 grupos taxonómicos. 
            
# Filtrado

En esta parte del programa se encuentra un loop que se encarga de filtrar los datos dependiendo del grupo taxonómico al que pertenecen **"phyllum"**, **"family"**, **"genus"**, o **"species"**. Esta sección del código tiene varios parámetros personalizables para poder satisfacer las necesidades de diferentes estudios. Uno de estos parámetros es el de la línea: 

```Rscript
#taxa_long_data[[i]]$Sample <- taxa_long_data[[i]]$Sample  %>% str_remove_all('[r\\d]') 
```
Esta línea permite agrupar los datos de las muestras para poner juntas las diferentes réplicas de cada muestra. Esta línea está comentada por default para evitar el agrupamiento de muestras con nombres similares. Ejemplo en un data frame con los siguientes datos: 

| Rank | TaxId | Scientific Name | SP1 | SP2 | SP3 | SP4
| --- | --- | --- | --- | --- | --- | --- |
| unclassified | 0 | Uknown | 68 | 1232 | 1696 | 22013 |
| superkingdom | 2 | Bacteria (bacteria) | 66708 | 76666 | 64937 | 25496 |
| genus | 131079 | Limnobacter| 0 | 0 | 0 | 90 |
| species | 2060312 | Altererythrobacter sp. B11 | 0 | 0 | 2 | 10 |

Al momento de hacer el filtrado de los datos no se van a reportar los resultados como muestras individuales "SP1, SP2, SP3, SP4", sino serán unidos todos los datos y serán reportados como las réplicas de SP. Es recomendable indicar códigos específicos para cada muestra y en caso de tener réplicas, poner los números. Es decir, en vez de poner: "Muestra1.1, Muestra1.2, Muestra2.1, Muestra2.2" , es mejor usar un código del estilo "A1, A2, B1, B2", así se agrupan las diferentes muestras con sus réplicas. En caso de no tener réplicas no es necesario hacer lo mencionado con anterioridad y dejar la línea comentada por default. Más ejemplos de cómo pueden ser los resultados pueden ser observados en [**Generación de Boxplots**](##Generación-de-Boxplots)


Dentro de este mismo loop se filtran los datos para poder obtener un **Top n** de los datos, es decir, obtener un Top n de especies, generos, familias y filos. Se puede modificar el tamaño del top dependiendo de las necesidades del estudio. El default de este programa es un Top 10. 

# Visualización de datos 

En esta sección se generan diferentes outputs en forma de documentos **.csv**, barplots y boxplots que ayudan a la visualización de los datos, así como son útiles para hacer análisis estadísticos posteriores debido a que ya pasaron por un normalizado y filtrado. 
## Generación de tablas 
## Generación de Boxplots 
