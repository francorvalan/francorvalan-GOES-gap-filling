# funciones_G2
# In[ ]:


#####--------------------------- ARCHIVO DE FUNCIONES GRUPO 2 Y 4  ####----------------------------


# In[1]:


#Importación de librerías
import matplotlib.pyplot as plt
import numpy as np

import rasterio
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.plot import show
from rasterio.fill import fillnodata
from rasterio.mask import mask

from scipy import signal


import geopandas as gpd
from geopandas.tools import sjoin

from shapely.geometry import Polygon , MultiPolygon, shape, Point # Lectura del shp y conversion de poligono 3D a 2D.

#import sys 
import os
import os.path

#%matplotlib notebook


# In[4]:


def convert_3D_2D(AOI):
    '''
    Esta función convierte Poligonos/Multipolígonos en 3D (Campo Z) a Polígonos/Multipolígonos 2D
    La función no presenta inconvenientes si el objeto de entrada se encuentre como objeto 2D.

    Parámetro:
        - AOI: ruta (path) del archivo vectorial del AOI
    Salida:
        - Lista con las geometrias del polígono en 2D.
    '''
    ROI = gpd.read_file(AOI)
    if ROI.crs.to_epsg() != 4326:
        ROI =  ROI.to_crs(4326)
    geometry = ROI.geometry
    
    
    new_geo = []
    for p in geometry:
        if p.has_z:
            if p.geom_type == 'Polygon':
                lines = [xy[:2] for xy in list(p.exterior.coords)]
                new_p = Polygon(lines)
                new_geo.append(new_p)
            elif p.geom_type == 'MultiPolygon':
                new_multi_p = []
                for ap in p:
                    lines = [xy[:2] for xy in list(ap.exterior.coords)]
                    new_p = Polygon(lines)
                    new_multi_p.append(new_p)
                new_geo.append(MultiPolygon(new_multi_p))
        
        else: new_geo = new_geo.append(p)
    return new_geo



# In[5]:


def Coregistration(infile, match, resampling=Resampling.nearest, return_objs = True, path='../data/proc/'):
    """
    Co-registro de Rasters con Rasterio (Alineación de celdas/grilla)

    Esta función reproyecta un archivo para que coincida con la forma y la proyección 
    de la grilla existente.
    La imagen corregistrada es guardada como .tif y tambien se retorna un array y sus metadatos..
    
    Parámetros de entrada: 
        - infile (str): ruta al archivo de entrada que se va a corregistrar.
        - match (str): ruta de la imagen de referencia con la forma y proyección deseadas.
        - resampling: metodo de remuestreo (ver a continuación)
        - return_objs: ¿Desea que devuelva un objeto tipo array?
        - path: directorio donde se guardaran las salidas

    Salida: 
        - dst_kwargs: metadatos del ráster corregistrado.
        - dst_array (array): Array del ráster corregistrado.
        - Escritura de archivo GeoTIF en ../proc/GFS_COR/"nombre_original" + "_COR"

    Metodologias de remuestreo/interpolación:
        - Resampling.nearest 
        - Resampling.bilinear
        - Resampling.cubic
        - Resampling.cubic_spline
        - Resampling.lanczos
        - Resampling.gauss
        - etc

    https://rasterio.readthedocs.io/en/latest/api/rasterio.enums.html#rasterio.enums.Resampling
    
    """
    outfile = path + "GFS_COR/" +'_COR'.join(os.path.splitext(os.path.basename(infile)))
    
    with rasterio.open(match) as reference:
        ref_meta = reference.meta

    # Lectura de la imagen a corregistrar
    with rasterio.open(infile) as src:
        dst_kwargs = src.meta.copy()
            
        # Configuración de los metaparámetros del salida (en funcion de la imagen de referencia)
        dst_kwargs.update({#"crs": ref_meta['crs'],
                           "crs": 'EPSG:4326', # Forzamos la asignación del EPSG ya que 
                           "transform": ref_meta['transform'],
                           "width": ref_meta['width'],
                           "height": ref_meta['height'],
                           "nodata": np.nan})



        # Escritura de la imagen a corregistrar
        with rasterio.open(outfile, "w+", **dst_kwargs) as dst:
            # En caso de presentar multibandas, itera para cada una de las bandas y escribe utilizando la función reproyectar
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=ref_meta['transform'],
                    dst_crs=ref_meta['crs'],
                    resampling=resampling)

            dst_array = dst.read(1) # Generación de un array de salida
        del dst   
    if return_objs:
        return dst_kwargs,dst_array
                


# In[6]:


def OutliersFilter(goes,umbral_outliers=230):
    ''' 
    Esta función aplica dos filtros para eliminar valores atípicos/erróneos de la imagen:

    - Filtro por umbral: eliminacióon de valores atipicos de temperatura, considerando aquellos que estén por debajo de valores históricos.
    - Filtro kernel: generación de un kernel 3x3 el cual retorna como salida la moda de los valores contíguos. Este filtro es aplicado 
                       ya que las imágenes goes presentan valores atípicos (muy bajos) dentro de la zona que presenta cobertura nubosa, pero
                       no son eliminados por por el filtro del umbral.

    Parámetros:
        - goes: datos de entrada, imagen a aplicar filtro.
        - umbral_outliers: temperatura mínima histórica.

    Salida: 
        - Imagen filtrada (array)

    '''
    # Filtro de valores atípicos de GOES por umbral (< ~ -40°C)
    goes[goes<umbral_outliers]=np.nan

    # Definición de tamaño de ventana para construcción de máscara
    kernel = np.zeros((3,3)) # genera valores 0 para datos con valores GOES y deja nan para datos vacios (nubes/filtro umbral)
    
    # Kernel de moda para generar una máscara, en caso de que la mayoría de la ventana sea nan, el resultado es nan.
    maskk = signal.convolve2d(goes, kernel,\
                              mode = 'same', boundary = 'fill', fillvalue = 0)

    # Finalmente la mascara (0: datos y NaN: datos) es interpolada para recuperar el borde de las áreas con datos.
    maskk=fillnodata(maskk, mask = (maskk == 0), max_search_distance = 3)

    goes_filt=maskk+goes
    
    return goes_filt


# In[7]:


# Funcion de rellenado de imagen GOES LST con GFS TSOIL
def GOES_fill(goes,gfs,gaps=True,outliers_filter=True,umbral_outliers=230, goes_nodata='nan', gfs_nodata=0 ):
    '''
    Esta función realiza un rellenado de datos vacios de la imagen GOES, utilizando los valores de GFS corregistrada.
    Para rellenar los datos, previamente se realiza un procesamiento de GFS, donde se rellenan las celdas vacias
    correspondientes a cuerpos de agua. Además se incoorpora la capacidad de eliminar datos atípicos/errores de GOES 
    con la función previamente definida (OutliersFilter).
    
    Post procesamiento de la imágen GOES y GFS
        - Eliminacion de outliers de GOES, setear la temperatura umbral en K
        - Relleno de huecos (gaps) de GFS
        - Relleno de GOES con GFS donde la primera no tiene datos
        - Generación de un producto de origen de datos: 1=GOES, 2=GFS, 0= sin datos en GOES & GFS
 
    Parámentros de Entrada:
        - gfs (array2D)
        - goes (array2D)
        - gaps (boolean): ¿Se desea rellenar los huecos de GFS?
        - outliers_filter (boolean): ¿Se desea eliminar valores atípicos/errores de GOES?
        - umbral_outliers= Temperatura Mínima en Kelvin para realizar el filtro por umbral.
        - goes_nodata: por omision 'nan', de lo contrario especificar
        - gfs_nodata: por omision 0 (cero), de lo cotrario especificar. Valores NaN expresar como 'nan'.

    Salida: 
        - orig_goes_rell (array): categorias de fuente de datos Nodata=0; GOES=1;GFS=2 
        - img_goes_rell (array): Imagen rellenada con valores de Temp en Kelvin.
    '''
    if goes_nodata != 'nan':
        goes[goes==goes_nodata] = np.nan

    if  gfs_nodata == 'nan':
        gfs[np.isnan(gfs)==1] = 0
    elif gfs_nodata != 0:
        gfs[gfs==gfs_nodata] = 0
    
    
    # Interpolación de datos nulos de GFS (cuerpos de agua) con IDW (Interpolacion Potencia Inversa a la Distancia)
    if gaps:
        fillnodata(gfs, mask = (gfs != 0), max_search_distance = 2)
        
    # Filtrado de valores atípicos por umbral y kernel
    if outliers_filter:
        goes = OutliersFilter(goes,umbral_outliers=230)
        
    # Rellenado de GOES=NaN con valores GFS
    img_goes_rell=goes.copy()
    
    img_goes_rell[(np.isnan(img_goes_rell)==1)&(gfs !=0)] = gfs[(np.isnan(img_goes_rell)==1)&(gfs !=0)]
    
    # Generación del raster de origen de pixeles
    orig_goes_rell=np.ones(goes.shape) # fuente 1 = GOES
    
    # valores de 2 para las celdas que contengan valores rellenados
    orig_goes_rell[np.isnan(goes)==1]=2 # fuente 2 = GFS
    
    # valores de NAN en celdas que no contengan valores en GOES ni GFS
    orig_goes_rell[(np.isnan(goes)==1) & (gfs==0)]=0 # fuente 0 = sin datos
    
    return orig_goes_rell,img_goes_rell


# In[9]:


def Post(AOI,filled_path, fill_source_path, Celsius=True): #new
    ''' 
    Esta función realiza un postprocesamiento de los datos:
        - Recorte de mascara para un AOI
        - Conversión de temperatura de Kelvin a Centígrados
        - Conversión de datos NAN a 9999
        - Redondeo de valores y conversion a enteros.

    Parámetros:
        - AOI (str): ruta (path) del archivo vectorial del AOI
        - filled_path (str): Directorio de la imagen GOES rellenada
        - fill_source_path (str): Directorio de la imagen SOURCE (fuente de origen de datos del relleno).
        - Celsius (boolean): ¿Desea convertir la temperatura a °Celsius?
    Salida: 
        - filled_masked (array): Array de la imagen GOES rellenada cortada con el aoi.
        - fill_source_masked (array): Array de la imagen SOURCE cortada con el aoi.
        - masked_transform: metadatos de transformación de las imágenes
    '''
#     if aoi.crs.to_epsg() != 4326:
#         aoi =  aoi.to_crs(4326)
    aoi = convert_3D_2D(AOI)
    
    filled = rasterio.open(filled_path, nodata='nan') #as out
    source = rasterio.open(fill_source_path, nodata=0) #as out

    #    out_array = out.read(0)
        
       
    filled_masked, masked_transform = mask(dataset=filled,shapes = aoi, crop = False,nodata=np.nan)    
    fill_source_masked, sour_mask_transf = mask(dataset=source,shapes = aoi, crop = False,nodata=0)        
    
    if Celsius:
        filled_masked=filled_masked-273.15 #int(round(273.15,0))     
        filled_masked[np.isnan(filled_masked)==1]=9999
        filled_masked=np.around(filled_masked).astype(np.int16)
    
    
    return filled_masked[0], fill_source_masked[0].astype(np.int8), masked_transform


# In[10]:



### Función para escribir los metadatos del producto #

def writeMetadata(GOES_path, GFS_path, fill_source_array, path='../data/out/'):
    '''
    Escribe los metadatos de los insumos utilizados y un resumen del resultado.
    Parametros de entrada:
    - GOES_path: string del directorio y nombre de la imagen GOES.
    - GFS_path: string del directorio y nombre de la imagen GFS.
    - fill_source_array: array con el producto de origen de datos (valores del array 0= sin datos, 1=GOES, 2=GFS)
    - path: ruta para almacenar el txt de los metadatos, por defecto 'path'.
    
    Salida: metadatos del producto LST GOES rellenados con TSOIL GFS.

     Estos metadatos aportan información sobre los insumos utilizados y brindan un resumen del resultado,siendo: 
    autor, 
        - fecha de la imagen GOES y GFS, 
        - porcentaje de píxeles GOES, 
        - porcentaje de píxeles rellenados con GFS, 
        - SRC y 
        - resolución espacial.
    '''

    
    
#    fill_source_array=fill_source_array[0]
    
    pixTot=fill_source_array.shape[0]*fill_source_array.shape[1]

    pixNAN=sum(sum(fill_source_array==0))
    pixVal=pixTot-pixNAN # pixeles con valores (sin espejos de agua)

    pOrig=f'Porcentaje de pixeles rellenados \
    {round(len(fill_source_array[fill_source_array==2])*100/pixVal,1)} %.'

    pRell=f'Porcentaje de pixeles originales \
    {round(len(fill_source_array[fill_source_array==1])*100/pixVal,1)} %.'

    f'Porcentaje de pixeles NAN {round(pixNAN*100/pixTot,1)} %.'
    
    GOES_GAD=GOES_path.split('_')[-3]

    GOES_GAD=GOES_path.split('_')[-1]

    if (GOES_GAD[0:2])=='GL':
        fechaGOES = f"Fecha Datos GOES: {GOES_GAD[8:10]}/\
{GOES_GAD[6:8]}/{GOES_GAD[2:6]} \
Hora: {GOES_GAD[10:12]}:00 (UTC-3)."
    else: fechaGOES = f"Fecha Datos GOES: Valor no hallado"
        
        
    GFS_GAD=GFS_path.split('_')[-1]        
    if (GFS_GAD[0:2])=='GL':        
        fechaGFS = f"Fecha Datos GFS: {GFS_GAD[8:10]}/\
{GFS_GAD[6:8]}/{GFS_GAD[2:6]} Hora Modelo: 00:00. \
Hora Pronostico: {GFS_GAD[10:12]}:00 (UTC-3)."
    else: fechaGFS = f"Fecha Datos GFS: Valor no hallado"


    archMeta=f'{GOES_path[-72:-4]}_META.txt' 
    file = open(path+ archMeta, "w")
    file.write("MAPA DE TEMPERATURAS EN SUPERFICIE."+ os.linesep)
    file.write("Datos LST GOES rellenados con  productos de pronóstico de temperatura a 2 m sobre el suelo de GFS, para las áreas con cobertura de nubes."+ os.linesep)
    file.write(''+ os.linesep)
    file.write("Autor: MAIE 2022-2023"+ os.linesep)
    file.write(""+ os.linesep)
    file.write('Valores de temperatura expresados en grados Celsius (ºC).'+ os.linesep)
    file.write(fechaGOES+ os.linesep)
    file.write(fechaGFS+ os.linesep)
    file.write(pOrig + os.linesep)
    file.write(pRell+ os.linesep)
    file.write('SRC: EPSG 4326.'+ os.linesep)
    file.write('Resolucion espacial: ~0.13º.'+ os.linesep)
    file.close()


# In[11]:


def write_GTiff(outfile, crs, transform, array, meta=None, nodata=None):
    '''
    Exportar imagen rellenada a formato Geotiff

    Parámetros: 
        - outfile (str): Directorio de salida para guardar el archivo 
        - crs (str EPSG):
        - transform (matriz afín):
        - array: array de imagen a guardar
        - meta: opcional, metadatos de la imagen a guardar
        - nodata: ¿que valor asignar a los nan al momento de escribir?
    Salida: 
        - Raster en formato GeoTIF
    '''
    if len(array.shape)==2:
        count=1
    else:
        count=array.shape[0]
    if not meta:
        meta = {}

    meta['driver'] = 'GTiff'
    meta['height'] = array.shape[-2]
    meta['width'] = array.shape[-1]
    meta['count'] = count
    meta['crs'] = crs
    meta['transform'] = transform

    if 'dtype' not in meta: #if no datatype is specified, use float32
        meta['dtype'] = np.float32
        
    if nodata==None:
        pass
    else:
        meta['nodata'] = nodata
    with rasterio.open(outfile, 'w', **meta) as dst:
        if count==1: #es una matriz bidimensional, la guardo
            dst.write(array.astype(meta['dtype']), 1)


# In[12]:


####### Grupo 4 severidad


def SevH(lst):
    """
Esta función categoriza los valores de temperatura en intensidad de helada:
  - mtr: matriz de temperatura en grados celcius sobre la cual trabaja la función, rasterio
Salida:
  - mtr: matriz categorizada, rasterio
"""
    mtr=lst.copy()
    for i in range(mtr.shape[0]):
        for j in range(mtr.shape[1]):
            
            if (mtr[i,j] >  0) & (mtr[i,j] != 9999) : # Sin helada
                mtr[i,j] = 1
            if mtr[i,j] <= 0 and mtr[i,j]> -2 : # Suaves
                mtr[i,j] = 2
            if mtr[i,j] <= -2 and mtr[i,j]> -4 :  # Moderadas
                mtr[i,j] = 3
            if mtr[i,j] <=-4 and mtr[i,j]>-6 : # Intensas
                mtr[i,j] = 4
            if mtr[i,j] <=-6 : # Muy intensas
                mtr[i,j] = 5

    return(mtr)


# In[ ]:





# In[ ]:



#################################################################
def ProcesarPares_Mirian(dir_root, nombreGOES, nombreGFS, verbose = True):

    # Path Harcodeados
    
    if verbose:
        print("*** Inicio de ProcesarPares_Mirian() ***")

    path_vectorial = f'{dir_root}data/raw_data/original_data/vectorial/'
    pathOUT = f'{dir_root}data/out/'
    pathOUT_sh = f'{dir_root}data/out/sh/'
    pathAOI = path_vectorial + 'aoi_shp.shp'
    ########################################################################
    
    #directorios de salida

    #para imagenes goes rellenadas 
    salida_goes_fill = pathOUT +'_FILL'.join(os.path.splitext(os.path.basename(nombreGOES)))
    salida_goes_source = pathOUT +'_SOURCE'.join(os.path.splitext(os.path.basename(nombreGOES)))

    #para imagenes goes rellenadas recortadas con area SANCOR
    nombrefill_mask=pathOUT +'_FILL_MASK'.join(os.path.splitext(os.path.basename(nombreGOES)))
    nombresource_mask = pathOUT +'_SOURCE_MASK'.join(os.path.splitext(os.path.basename(nombreGOES)))

    # para salida de severidad de heladas
    nombreSH=pathOUT_sh+'SH_'+os.path.basename(nombreGOES)[58:68]+'.tif' #+'_'+os.path.basename(nombreGOES)[66:68]+'.tif'

    if verbose:
        print("Open GOES: ", nombreGOES)

    # Carga de imagen GOES recortada y guardado de metadatos para luego exportar.
    GOES_r= rasterio.open(nombreGOES)
    gt_GOES_r = GOES_r.transform # geotransformación.
    crs_GOES_r   = GOES_r.crs  # crs de la capa.
    GOES_r_band= GOES_r.read(1) # Array de la banda 1.

        # 1. CORREGGISTRO
    if verbose:
        print("Open GFS:", nombreGFS)
        print("Inicio de Corregistro")

    gfs_kwargs,gfs_array=Coregistration(infile = nombreGFS, 
                 match= nombreGOES,
                 resampling=Resampling.nearest
               )

    # 2. RELLENO DE DATOS 
    if verbose:
        print("Inicio de Rellenado")
    source,GOES_filled=GOES_fill(GOES_r_band,gfs_array)


    # 3. EXPORTACION DE IMAGEN GOES RELLENADA Y FUENTE DE DATOS

    if verbose:
        print("Carpeta OUT:", pathOUT)

    if verbose:
        print("SAVE 1 de 6 - Imagen GOES rellenada")
    write_GTiff(salida_goes_fill,crs_GOES_r,gt_GOES_r,GOES_filled,nodata=9999)

    if verbose:
        print("SAVE 2 de 6 - Fuente de datos")
    write_GTiff(salida_goes_source,crs_GOES_r,gt_GOES_r,source,nodata=0)

    # 4. POST PROCESO: CORTE DE AREA DE INTERES SANCOR, CONVERSION A TEMPERATURA CELSIUS Y REDONDEO DE VALORES DE TEMPERATURA
    if verbose:
        print("   Corte del area de interes")
    fill_mask, fill_source_mask, mask_transform=Post(pathAOI,salida_goes_fill, salida_goes_source)

    # 5. CREACION DE METADATO DE IMAGEN RELLENADA 
    #Estos metadatos aportan información sobre los insumos utilizados y brindan un resumen del resultado, siendo: autor, fecha de la imagen GOES y GFS, porcentaje de píxeles GOES, porcentaje de píxeles rellenados con GFS, SRC y resolución espacial.

    if verbose:
        print("SAVE 3 de 6 - Metadato de la Goes Rellenado")
    writeMetadata(nombreGOES, nombreGFS, fill_source_mask)


    # 6. CLASIFICACION DE SEVERIDAD DE HELADA
    #categorización de los valores de temperatura en intensidad de helada
    if verbose:
        print("   Clasificacion de Severidad de helada")
    severidad=SevH(fill_mask)


    # Exportación de salidas

    if verbose:
        print("Save 4 de 6 - Helada Recortado")
    write_GTiff(nombreSH,crs_GOES_r,mask_transform,severidad,nodata=9999) # Helada Recortado

    if verbose:
        print("Save 5 de 6 - GOES rellenado cortado")
    #write_GTiff(nombrefill_mask,crs_GOES_r,mask_transform,source,nodata=9999) # GOES rellenado cortado
    write_GTiff(nombrefill_mask,crs_GOES_r,mask_transform,fill_mask,nodata=9999) # GOES rellenado cortado

    if verbose:
        print("Save 6 de 6 - Fuente de Datos Recortado")
    write_GTiff(nombresource_mask,crs_GOES_r,mask_transform,fill_source_mask,nodata=9999) # Fuente de datos (goes o gfs pixel a pixel)
    

    if verbose:
        print("*** Fin de ProcesarPares_Mirian() ***")


    return


############## Función de activación ###################

def proceso_grupo2(dir_root = '../', verbose=True):
    import funciones_G1
    w_goes = f'{dir_root}data/raw_data/work_data/GOES/'
    w_gfs = f'{dir_root}data/raw_data/work_data/GFS/'
    path_out = f'{dir_root}data/out/'

    list_1, list_2 = funciones_G1.lista_procesar_proc(path_tif1=w_goes,path_tif2=w_gfs,path_final=path_out)
    if list_1 != []:
        for file_a, file_b in zip(list_1, list_2):
            goes = w_goes+file_a
            gfs = w_gfs+file_b
            ProcesarPares_Mirian(dir_root, goes, gfs, verbose = True)
    else:
        if verbose==True:
            print('No hay imagenes para procesar')
    return















