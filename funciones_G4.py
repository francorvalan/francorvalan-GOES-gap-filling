import rasterio
import pandas as pd
import os
import datetime
import numpy as np


def matcheo_hha(dir_root = '../', folder = 'sh'):
    '''
    Crea una lista de imagenes con sus fechas de las carpetas sh y hha. Sirve para la función control_hha.
    folder: 'sh' para buscar en la carpeta sh. 'hha' para buscar en la carpeta hha
    '''
    ext = '.tif'
    ubi_date = -1
    prefix = 0
    suffix = 10
    format='%Y%m%d%H'

    if folder == 'sh':
        path = f'{dir_root}data/out/sh/'
        prejijo = 'SH_'
    elif folder == 'hha':
        path = f'{dir_root}data/out/hha/'
        prejijo = 'HHA_'

    path_files = [f for f in os.listdir(path) if f.endswith(ext)]
    path_files = [x for x in path_files if prejijo in x]

    # Itera en todos los nombres de la carpeta original
    list_date = []
    for file in path_files:
        # Toma el nombre del archivo .nc original del GOES
        fecha = file.split("_")[ubi_date][prefix:suffix] 
        fecha_utc = pd.to_datetime(fecha, format=format, utc=True) # se interpreta en formato de fecha en UTC
        # Guarda la fecha UTC en str en formato juliano y gregoriano
        utc = f'{fecha_utc:%Y%m%d%H}' 
        # Guarda en str la fecha
        img_ogoes = f'{utc}'
        list_date.append(img_ogoes)

    df_imgs = pd.DataFrame({'fecha_img': list_date, 'nom_archivo': path_files})
    
    return df_imgs


def control_hha(verbose = True):
    '''
    Controla disponibilidad de imagenes para procesar.
    Genera una lista de horas para calcular las heladas del último día y 
    compara con las imágenes ubicadas en la carpeta sh. 
    Si faltan, avisa. Sino, entrega la lista de imagenes a procesar. 
    '''
       
    start = datetime.datetime.utcnow() - datetime.timedelta(days=1) ### Cambiar
    end = datetime.datetime.utcnow()

    start = f'{start:%d-%m-%Y} 12'
    end = f'{end:%d-%m-%Y} 11'

    DATES = pd.date_range(start, end, freq="H")

    heladas_en_carpeta = matcheo_hha(folder='sh') #lista_matcheo(path = path_hha, source='prc')

    horas_acumuladas = []
    for DATE in DATES:
        d = f'{DATE:%Y%m%d%H}'
        horas_acumuladas.append(d)

    df_horas_acumuladas = pd.DataFrame({'fecha_acum':horas_acumuladas})

    procesar = pd.merge(df_horas_acumuladas, heladas_en_carpeta[["fecha_img", 'nom_archivo']], left_on='fecha_acum', right_on='fecha_img', how='left')
    procesar_null = procesar[pd.isnull(procesar.nom_archivo)]
    procesar_notull = procesar[pd.notnull(procesar.nom_archivo)]
    list_procesar = list(procesar_notull.nom_archivo)

    if len(pd.isnull(list(procesar_null.nom_archivo))*1) > 0:
        if verbose == True:
            print(f'[!] Faltan las imagenes de los días {list(procesar_null.fecha_acum)} para poder realizar la tarea')
        list_procesar = []
        
    return list_procesar


def guardar_GTiff(fn, crs, transform, mat, meta=None, nodata=None, bandnames=[]):
    if len(mat.shape)==2:
        count=1
    else:
        count=mat.shape[0]

    if not meta:
        meta = {}

    meta['driver'] = 'GTiff'
    meta['height'] = mat.shape[-2]
    meta['width'] = mat.shape[-1]
    meta['count'] = count
    meta['crs'] = crs
    meta['transform'] = transform

    if 'dtype' not in meta: #if no datatype is specified, use float32
        meta['dtype'] = np.float32
    

    if nodata==None:
        pass
    else:
        meta['nodata'] = nodata

    with rasterio.open(fn, 'w', **meta) as dst:
        if count==1: #es una matriz bidimensional, la guardo
            dst.write(mat.astype(meta['dtype']), 1)
            if bandnames:
                dst.set_band_description(1, bandnames[0])
        else: #es una matriz tridimensional, guardo cada banda
            for b in range(count):
                dst.write(mat[b].astype(meta['dtype']), b+1)
            for b,bandname in enumerate(bandnames):
                dst.set_band_description(b+1, bandname)#   


def hha(dir_root='../', verbose=True):
    """
    Esta función toma el producto de severidad de helada desde las 12:00 del dia anterior a las 11:00 del presente dia y
    devuelve un tiff que indica en cuales pixels se presentaron heladas acumuladas:
    - path: dirección donde se encuentran los productos de severidad de helada, string
    - driver: formato de salida del producto, string
    Salida:
    - map: matriz calculada, rasterio
    """
    
    path_sh = f'{dir_root}data/out/sh/'
    path_hha = f'{dir_root}data/out/hha/'
    
    img_en_carpeta = True

    # Control de no repetición de la tarea
    if list(matcheo_hha(folder='hha').fecha_img.values) != []: # no existen tantas imagenes como fechas de heladas para procesar
        if control_hha(verbose=False) != []:
            img_en_carpeta = True
            if sum(matcheo_hha(folder='hha').fecha_img.values == control_hha()[0].split('_')[-1][0:10]*1) > 0: # si ya estaán creadas las hha:
                img_en_carpeta = False
    if control_hha() == [] :
        img_en_carpeta = False

    if img_en_carpeta:
        bands=[]
        
        lista_hha = control_hha(verbose=False) # Recibe la lista de imagenes sólo si el conjunto está completo
        
        for img_name in lista_hha:
            img= rasterio.open(path_sh + img_name)
            img_bands = img.read()
            img_1=img_bands[0]

            img_b = img_1.copy()
            img_b[img_1 > 1] = 1
            img_b[img_1 == 1] = 0
            img_b[img_1 == 9999] = np.nan

            # for x in range(img_1.shape[0]):
            #     for y in range(img_1.shape[1]):
            #         if img_1[x,y] > 1: #Cuidado este valor es simulado , cambiar a 1 (lo cambié. german)
            #             img_1[x,y] = 1
            #         else:
            #             img_1[x,y] = 0

            bands.append(img_b)
            
        map=(bands[0]+bands[1]+bands[2]+bands[3]+bands[4]+bands[5]+bands[6]
            +bands[7]+bands[8]+bands[9]+bands[10]+bands[11]+bands[12]
            +bands[13]+bands[14]+bands[15]+bands[16]+bands[17]+bands[18]
            +bands[19]+bands[20]+bands[21]+bands[22]+bands[23])
        
        map[map == np.nan] = 9999
        # Se define una variable para almacenar la información del metadato del objeto Rasterio
        meta=img.meta

        # Se actualiza(sobreescribe) el objeto meta, para definir el tipo de dato como entero,nodata=0
        # Se genera una varaible que almacena  el CRS en el cual se esta trabajando
        # Se genera una varaible con que almacena el transform del objeto Rasterio sobre el cual se esta trabajando

        # Se genera el nombre del archivo a gardar
        tif_name = 'HHA_'+lista_hha[0].split('_')[-1][0:10]+'-'+lista_hha[-1].split('_')[-1][0:10]+'.tif'

        # Se define una variable para almacenar la información del metadato del objeto Rasterio
        meta=img.meta

        # Se actualiza(sobreescribe) el objeto meta, para definir el tipo de dato como entero,nodata=0
        meta.update(dtype=rasterio.uint8,width=img_1.shape[0] ,height=img_1.shape[1],nodata=9999,count= 1)

        # Se genera una varaible que almacena  el CRS en el cual se esta trabajando
        crs=img.crs

        # Se genera una varaible con que almacena el transform del objeto Rasterio sobre el cual se esta trabajando
        transform=img.transform

        # Se exporta la matriz recategorizada como archivo Tiff
        guardar_GTiff(path_hha+tif_name, crs, transform, map)

        if verbose == True:
            print(f'Procesando img: {tif_name}')
    else:
        if verbose == True:
            print('[!] No hay imagenes nuevas para crear')