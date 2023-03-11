# TFM Bioestadística y Bioinformática

Trabajo de final de máster de Bioinformática y Bioestadística para la Universitat Oberta de Catalunya. Aquí se hace un pequeño resumen de las etapas a realizar durante este trabajo. El objetivo principal de este trabajo es describir un *workflow* que permita la descripción y lanzamiento de tareas que analicen cualquier conjunto de datos ómicos de diferente índole (metabolómicos, transcriptómicos...). Para ello se puede dividir este proyecto en las siguientes etapas:

1. Hacer una revisión de las aplicaciones similares al objetivo del TFM
2. Revisar los *pipelines* en la literatura para datos ómicos de diferente tipo y que librerías usar
3. *Quality Control* (QC) y normalización de de los datos
4. Escribir scripts que hagan análisis descriptivos de los datos específicos para cada tipo de dato ómico
5. Escribir scripts de *feature selection* y *Machine Learning*
6. Integrar esto en una aplicación escrita con librería *Shiny*.

## 1. Revisión de aplicaciones disponibles

### 1.1. [Metaboanalyst 5.0](https://www.metaboanalyst.ca/)

Esta aplicación web esta especializada en el análisis de datos metabolómicos desde *raw data* hasta inferencia estadística. Es de las aplicaciones más usadas en metabolómica

### 1.2. [SIMCA](https://www.sartorius.com/en/products/process-analytical-technology/data-analytics-software/mvda-software/simca)

### 1.3. [XCMS](https://xcmsonline.scripps.edu/landing_page.php?pgcontent=mainPage)

### 1.4. [mixOmics](http://mixomics.org/)

Este paquete de R tiene diversas funcionalidades para integrar datos multiómicos, ya sean diferentes muestras (*batches*) o diferentes análisis de las mismas muestras (RNAseq, metabolómicos...) para integrarlos y aplicar algoritmos ML supervisados y no supervisados. Esta libería entonces tiene objetivos muy parecidos a este proyecto. Sin embargo, mixOmics espera una serie de datos normalizados. 

### 1.5. [Galaxy](https://usegalaxy.org/)

### 1.6. [Genevestigator](https://genevestigator.com/)


