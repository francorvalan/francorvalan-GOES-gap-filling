###########
# V1.0
# Unificación 15/12/2022
###########

# Librerias
import os
import s3fs
from rptree.rptree import DirectoryTree
import rioxarray as rxr
import calendar
import datetime
import re
from dateutil import tz
from pyproj import CRS
import xarray as xr
import pandas as pd
from os.path import isfile, join
import gc
import matplotlib.pyplot as plt
import numpy as np
import copy
import sys
import requests
import warnings
import time

import schedule


# time.sleep()

# Función de seteo de directorio raiz y carpetas

def crear_carpetas(dir_root = '../'):
    # Make new directories to store files
    
    if not os.path.isdir(dir_root + "data"):
        os.mkdir(dir_root + "data")
    if not os.path.isdir(dir_root + "data/" +"out"): 
        os.mkdir(dir_root + "data/" +"out")
    if not os.path.isdir(dir_root + "data/" +"out/" + 'sh'): 
        os.mkdir(dir_root + "data/" +"out/" + 'sh')
    if not os.path.isdir(dir_root + "data/" +"out/" + 'hha'): 
        os.mkdir(dir_root + "data/" +"out/" + 'hha')
    if not os.path.isdir(dir_root + "data/" +"proc"):
        os.mkdir(dir_root + "data/" +"proc")
    if not os.path.isdir(dir_root + "data/" +"proc/" + 'GFS_COR'):
        os.mkdir(dir_root + "data/" +"proc/" + 'GFS_COR')
    if not os.path.isdir(dir_root + "data/" +"raw_data"):
        os.mkdir(dir_root + "data/" + "raw_data")
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "original_data"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "original_data")
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "original_data/" + "GFS"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "original_data/" + "GFS")
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "original_data/" + "GOES"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "original_data/" + "GOES")
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "original_data/" + "vectorial"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "original_data/" + "vectorial")
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "work_data"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "work_data")
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "work_data/" + "GFS"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "work_data/" + "GFS")    
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "work_data/" + "GOES"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "work_data/" + "GOES")
    if not os.path.isdir(dir_root + "data/" +"raw_data/" + "work_data/" + "vectorial"):
        os.mkdir(dir_root + "data/" +"raw_data/" + "work_data/" + "vectorial")

    return

    # Función de nombrado de archivo de salida GOES

# # # # # # # # # Funciones relacionadas a fechas...
def JulianDate_to_GregorianDate(julian_date):
    # https://stackoverflow.com/questions/13943062/extract-day-of-year-and-julian-day-from-a-string-date
    
    #############################
    y = int(julian_date[0:4])
    jd = int(julian_date[4:7])
    the_hour = julian_date[7:9]
    #############################
    
    month = 1
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month += 1

    date = datetime.date(y,month,jd).strftime("%Y/%m/%d")

    # Special Selection
    date = re.sub("/","",date)
    date = date + the_hour

    return date

def GUTC_to_GLocal(gregorian_date):
    
    final = ":00:00Z"
    nuevo_formato = gregorian_date[0:4] + "-" + gregorian_date[4:6] + "-" + gregorian_date[6:8] + "T" + gregorian_date[8:10] + final
    from_zone = tz.gettz('UTC')
    to_zone = tz.gettz('America/Argentina/Buenos_Aires')
    # json_data = {'time': "2021-10-08T08:17:42Z"} 
    json_data = {'time': nuevo_formato} 
    # json_data['time']
    utc = datetime.datetime.strptime(json_data['time'], "%Y-%m-%dT%H:%M:%SZ")
    utc = utc.replace(tzinfo=from_zone)
    local = utc.astimezone(to_zone)
    local = str(local)

    # Special Sentences
    local = re.sub("-","",local)
    local = re.sub(" ","",local)
    local = local[0:10]

    
    return(local)


def GLocal_to_JLocal(glocal):
    recambio = glocal[0:4] + "-" + glocal[4:6] + "-" + glocal[6:8]
    jlocal = datetime.datetime.strptime(recambio, "%Y-%m-%d").strftime('%Y%j')
    jlocal = jlocal + glocal[8:10]
    return jlocal

def FullDateGenerator_from_JD(julian_download_date):
    goes_julian_download_date = julian_download_date

    goes_gregorian_download_date = JulianDate_to_GregorianDate(julian_date = goes_julian_download_date)

    glocal = GUTC_to_GLocal(gregorian_date = goes_gregorian_download_date)

    jlocal = GLocal_to_JLocal(glocal)

    full_date = "JAD" + goes_julian_download_date + "_GAD" + goes_gregorian_download_date + "_JL" + jlocal + "_GL" + glocal
    return full_date


def GFecha_hoyUTC():
    fecha_hoy = datetime.datetime.utcnow() # Fecha y Hora UTC para descargar!
    fecha_hoy = str(fecha_hoy)
    fecha_hoy = fecha_hoy[0:4] + "-" + fecha_hoy[5:7] + "-" + fecha_hoy[8:10]
    return fecha_hoy

def GFecha_ayerUTC():
    fecha_hoy = datetime.datetime.utcnow()
    ayer = fecha_hoy - datetime.timedelta(days = 1)
    str_fecha_hoy = str(fecha_hoy)
    str_ayer = str(ayer)[0:10]
    return str_ayer

def GFechaHora_descarga_utc_hoy_GOES():
    las_horas = [f'{number:02}' for number in range(24)]
    fecha_hoy = GFecha_hoyUTC()
    gd_hora_descarga_goes = [fecha_hoy + " " + x for x in las_horas]
    return gd_hora_descarga_goes

def GFechaHora_descarga_utc_ayer_GOES():
    las_horas = [f'{number:02}' for number in range(24)]
    fecha_ayer = GFecha_ayerUTC()
    gd_hora_descarga_goes = [fecha_ayer + " " + x for x in las_horas]
    return gd_hora_descarga_goes

def GFechaHora_Ingresa_Usuario_Descargar_GFS():
    fecha_hoy = GFecha_hoyUTC()
    return fecha_hoy



def Ayer_UTC_GFS():
    fecha_hoy = datetime.datetime.utcnow()
    ayer = fecha_hoy - timedelta(days = 1)
    return fecha_hoy

def GFechaHora_Archivos_para_DescargarGOES():
    fechas = GFechaHora_Ingresa_Usuario_Descargar_GOES()
    fechas = [x.replace(' ', '').replace('-', '') for x in fechas]
    return fechas

def GFechaHora_Archivos_Server_GOES(dir_root = "../"):
    control_folder_goes = dir_root + "data/raw_data/original_data/GOES/"
    archivos_original_goes = os.listdir(control_folder_goes)
    fechas_archivos_goes = [x[23:32] for x in archivos_original_goes]
    return fechas_archivos_goes



# # # # # # # # # # # # # # # # # # # Funciones Intermedias
def sufijo_tif(file, source):
    
    if source == 'goes':
        ubi_date = -3
        prefix = 1
        suffix = 10
        format='%Y%j%H'
    elif source == 'gfs':
        ubi_date = -1
        prefix = 0
        suffix = 10
        format='%Y%m%d%H'
    elif source == 'prc':
        ubi_date = -3
        prefix = 3
        suffix = 13
        format='%Y%m%d%H'

    # Toma el nombre del archivo .grib original del GFS
    fecha_jul = file.split("_")[ubi_date][prefix:suffix] 
    fecha_utc = pd.to_datetime(fecha_jul, format=format, utc=True) # se interpreta en formato de fecha en UTC

    # Guarda la fecha UTC en str en formato juliano y gregoriano
    j_utc = f'{fecha_utc:JAD%Y%j%H}' 
    g_utc = f'{fecha_utc:GAD%Y%m%d%H}' 

    # Convierte la fecha UTC a hora local
    fecha_local = fecha_utc.tz_convert('America/Cordoba')

    # Guarda la fecha local en str en formato juliano y gregoriano
    j_local = f'{fecha_local:JL%Y%j%H}'
    g_local = f'{fecha_local:GL%Y%m%d%H}'

    # Guarda en str el nombre completo del archivo tiff
    file_name = f'{j_utc}_{g_utc}_{j_local}_{g_local}.tif'

    return file_name


def lista_procesar(original_path, work_path, source):

    if source == 'goes':
        ext = '.nc'
        ubi_date = -3
        prefix = 1
        suffix = 10
        format='%Y%j%H'
    elif source == 'gfs':
        ext = '.grib'
        ubi_date = -1
        prefix = 0
        suffix = 10
        format='%Y%m%d%H'
    elif source == 'prc':
        ext = '.tif'
        ubi_date = -3
        prefix = 3
        suffix = 13
        format='%Y%m%d%H'

    original_data_path = [f for f in os.listdir(original_path) if f.endswith(ext)]

    # Itera en todos los nombres de la carpeta original
    list_original_date = []
    for file in original_data_path:
        # Toma el nombre del archivo .nc original del GOES
        fecha_ogoes = file.split("_")[ubi_date][prefix:suffix] 
        fecha_ogoes_utc = pd.to_datetime(fecha_ogoes, format=format, utc=True) # se interpreta en formato de fecha en UTC
        # Guarda la fecha UTC en str en formato juliano y gregoriano
        og_utc = f'{fecha_ogoes_utc:%Y%m%d%H}' 
        # Guarda en str la fecha
        img_ogoes = f'{og_utc}'
        list_original_date.append(img_ogoes)


    work_data_path = [f for f in os.listdir(work_path) if f.endswith('.tif')]

    # Itera en todos los nombres de la carpeta destino
    list_work_date = []
    for file in work_data_path:
        # Toma el nombre del archivo .tif procesado de GOES
        fecha_wgoes = file.split("_")[-3][3:]
        fecha_wgoes_utc = pd.to_datetime(fecha_wgoes, format='%Y%m%d%H', utc=True) # se interpreta en formato de fecha en UTC
        # Guarda la fecha UTC 
        wg_utc = f'{fecha_wgoes_utc:%Y%m%d%H}' 
        # Guarda en str ela fecha
        img_wgoes = f'{wg_utc}'
        list_work_date.append(img_wgoes)

    d1 = pd.DataFrame({'fecha_img': list_original_date, 'nom_archivo_ori': original_data_path})
    d2 = pd.DataFrame({'fecha_img': list_work_date, 'nom_archivo_prc': work_data_path})
    procesar = pd.merge(d1, d2[["nom_archivo_prc",'fecha_img']], left_on='fecha_img', right_on='fecha_img', how='left')
    procesar = procesar[pd.isna(procesar.nom_archivo_prc)]
    procesar = list(procesar.nom_archivo_ori)
    return procesar


def Detectar_Procesamiento(input_files, output_files, pos_input, pos_output):
        # Listamos los archivos en input y en output.
        # input_files = os.listdir(input_folder)
        # output_files = os.listdir(output_folder)

        # Del input, tomamos una porcion del nombre
        # pos_input = [13, 23]
        gfd_input = [x[pos_input[0]:pos_input[1]] for x in input_files]

        # Del output, tomamos una porcion del nombre
        # pos_output = [30, 40]
        gfd_output = [x[pos_output[0]:pos_output[1]] for x in output_files]

        # Determinamos si cada input esta en la carpeta output.
        # Si no esta en la carpeta otorga un True, lo que significa
        # que el preprocesamiento debe ser ejecutado sobre esa imagen.
        dt_procesar = [x not in gfd_output for x in gfd_input]

        return dt_procesar
    
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # GOES
def NombrarGOES(original_file_name):
    # original_file
    # goes_julian_download_date = original_file_name[23:32]
    julian_download_date = original_file_name[23:32]

    full_date = FullDateGenerator_from_JD(julian_download_date)

    new_file_name = "GOES_" + original_file_name[3:14] + "_" + full_date + ".tif"

    return(new_file_name)

############# Descarga de GOES ###################
def descarga_goes(start, end, verbose=True, dir_root = '../'):
    
    # Seteo de direcciones de apertura y guardado
    goes_ofolder = f'{dir_root}data/raw_data/original_data/GOES/'

    # Seteo satélite y producto a descargar de GOES
    bucket_name = 'noaa-goes16' # Cambiar por 'noaa-goes17' para GOES-17
    product_name = 'ABI-L2-LSTF'# Otros productos GOES-R ABI L2+ ('https://docs.opendata.aws/noaa-goes16/cics-readme.html'

    lista_descarga_goes = control_descarga(start=start,end=end,path=goes_ofolder, source='goes')

    if lista_descarga_goes != []:
        
        # Inicio de sesión anónimo en AWS
        fs = s3fs.S3FileSystem(anon=True)
            
        if verbose == True:
            print('** Iniciando la descarga de GOES **')

        # Se genera una seria por hora entre las fechas seleccionadas, que todavía no han sido descargadas
        DATES = pd.to_datetime(lista_descarga_goes, format= '%Y%m%d%H')

        # Itera en todas las fechas
        for DATE in DATES:
            # Toma una fecha y busca el archivo .nc que conicida en la base de AWS
            files = fs.ls(f"{bucket_name}/{product_name}/{DATE:%Y/%j/%H/}", refresh=True)          

            if files != []:
                # Toma la última en el caso de que haya más de una imagen de la fecha
                last_file = files[0]

                # Se crea la dirección de guardado, el .nc se guarda con el nombre original
                path = goes_ofolder+last_file.split("/")[-1]
                
                # Descarga
                fs.download(last_file, path)
                
                if verbose == True:
                    print(f'    - Imagen del {DATE:%Y-%m-%d %H:00} UTC+0 descargada correctamente.')
            else:
                if verbose == True:
                    print(f'    - [!] La imagen del {DATE:%Y-%m-%d %H:00} UTC+0 no está disponible para la descarga.')
    else:
        if verbose == True:
            print('[!] No hay imagenes nuevas de GOES para descargar')    
    return


def procesamiento_goes(dir_root = '../', celsius = False, verbose=True):
    # Seteo de direcciones de apertura y guardado
    goes_ofolder = f'{dir_root}data/raw_data/original_data/GOES/'
    goes_wfolder = f'{dir_root}data/raw_data/work_data/GOES/'

    if verbose == True:
        print(f'** Inicio de procesamiento GOES **')
        if lista_procesar(goes_ofolder, goes_wfolder, source='goes') == []:
            print('  - [!] No hay nuevas imagenes para procesar')

    for file in lista_procesar(goes_ofolder, goes_wfolder, source='goes'):

        # Obtiene la dirección del archivo .nc
        pathfiles = f'{goes_ofolder}{file}'
        
        # Abro el .nc
        nc = rxr.open_rasterio(pathfiles, mask_and_scale=True)

        # Se establece el la referencia espacial y se define en el archivo
        spatial_ref = 'PROJCS["unnamed",GEOGCS["unknown",DATUM["unnamed",SPHEROID["Spheroid",6378137,298.2572221]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]]],PROJECTION["Geostationary_Satellite"],PARAMETER["central_meridian",-75],PARAMETER["satellite_height",35786023],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],EXTENSION["PROJ4","+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x"]]'
        nc_ref = nc.rio.set_crs(spatial_ref)

        nc_ref = nc_ref.rio.reproject("EPSG:4326")

        # Clip por bbox
        nc_clip = nc_ref.rio.clip_box(
            minx=-68,
            miny=-42,
            maxx=-55,
            maxy=-21,
            crs="EPSG:4326",
        ).astype(float)

        nc_clip.rio.write_crs("EPSG:4326", inplace=True)

        # Cambio de Kelvin a Grados Celsius
        if celsius == True:
            nc_clip = nc_clip - 273.15

        # Se crea el nombre del archivo y dirección del archivo a guardar
        file_name = sufijo_tif(file, source='goes')
        tiff_name = f'{goes_wfolder}GOES_ABI-L2-LSTF_{file_name}'

        # Se guarda el .tif
        nc_clip.LST.rio.to_raster(tiff_name)
    
        if verbose == True:
            ffile=file_name.split('_')[-3][3:]
            fecha_str = pd.to_datetime(ffile, format='%Y%m%d%H')
            print(f'    - Imagen {fecha_str:%Y-%m-%d %H:00} UTC+0 procesada correctamente.')
    
    return


def Detect_Descarga_Archivos_GOES(gdh_usuario_descarga_goes, gdh_en_carpeta01, dir_root = "../"):
    gdh_en_descarga = [x.replace("-", "").replace(" ", "") for x in gdh_usuario_descarga_goes]
    jdh_en_descarga = [GLocal_to_JLocal(x) for x in gdh_en_descarga]

    dt_descarga01 = [x in gdh_en_carpeta01 for x in gdh_en_descarga]
    dt_descarga01 = [x == False for x in dt_descarga01]
    return dt_descarga01


def DescargarGoesMOD(start, verbose = True, dir_root = "../"):
    # Seteo satélite y producto a descargar de GOES
    bucket_name = 'noaa-goes16' # Cambiar por 'noaa-goes17' para GOES-17
    product_name = 'ABI-L2-LSTF'# Otros productos GOES-R ABI L2+ ('https://docs.opendata.aws/noaa-goes16/cics-readme.html'

    # Inicio de sesión anónimo en AWS
    fs = s3fs.S3FileSystem(anon=True)

    # Julian Date Format - Original
    jdf_original = start
    jdf_new = jdf_original[0:4] + "/" + jdf_original[4:7] + "/" + jdf_original[7:10]
    
    gdf_original = JulianDate_to_GregorianDate(jdf_original)
    gdf_new = gdf_original[0:4] + "/" + gdf_original[4:6] + "/" + gdf_original[6:8] + " " + gdf_original[8:10]
    
    full_date = FullDateGenerator_from_JD(julian_download_date = jdf_original)
    armado = jdf_new + " *** " + gdf_new
    
    # fecha_formato_nuevo
    busqueda_especial01 = bucket_name + "/" + product_name + "/" +  jdf_new

    files = fs.ls(busqueda_especial01, refresh=True)

    # Si hay al menos una imagen para esa fecha y hora...
    if len(files) > 0: 
        last_file = files[0]

        nombre_original_goes = last_file.split("/")[-1]
        # nuevo_nombre = nombre_original[0:22] + full_date + ".nc"
        # Se crea la dirección de guardado, el .nc se guarda con el nombre original
        path = dir_root + "data/raw_data/original_data/GOES/" + nombre_original_goes
        # path = dir_root + "data/raw_data/original_data/GOES/" + nuevo_nombre
        path

        # Descarga
        fs.download(last_file, path)
        
        
        if verbose == True:
             print(f'    - Descarga Img GOES *** UTC +00 *** {armado} - Descargada Correctamente.')

    else:
        if verbose == True:
            print(f'    - Descarga Img GOES *** UTC +00 *** {armado} - Sin producto en el servidor.')

            
def Descarga_Controlada_GOES(dir_root = "../", gdh_usuario_descarga_goes = GFecha_hoyUTC(), verbose = True):
    
    dir_root = "../"
    control_folder_goes = dir_root + "/" + "data/raw_data/original_data/GOES/"
    archivos_original_goes = os.listdir(control_folder_goes)

    jdh_en_carpeta01 = os.listdir(control_folder_goes)
    jdh_en_carpeta01 = [x[23:32] for x in jdh_en_carpeta01]
    
    # gdh_usuario_descarga_goes = GFechaHora_Ingresa_Usuario_Descargar_GOES()
    gdh_en_descarga = [x.replace("-", "").replace(" ", "") for x in gdh_usuario_descarga_goes]
    jdh_en_descarga = [GLocal_to_JLocal(x) for x in gdh_en_descarga]
    
    # Vamos a comparar el juliano de descarga con los julianos de nuestra carpeta
    # Los julianos de descarga que no estan en la carpeta, se los debe intentar descargar.
    dt_descarga01 = [x in jdh_en_carpeta01 for x in jdh_en_descarga]
    dt_descarga01 = [x == False for x in dt_descarga01]
    
    for x in range(len(dt_descarga01)):
        if dt_descarga01[x]: 
            DescargarGoesMOD(start = jdh_en_descarga[x], verbose = verbose, dir_root = "../")
        else: 
            fecha_formato_nuevo_juliano = jdh_en_descarga[x][0:4] + "/" + jdh_en_descarga[x][4:7] + "/" + jdh_en_descarga[x][7:10]
            armado = fecha_formato_nuevo_juliano + " *** " + gdh_usuario_descarga_goes[x] 
            
            if verbose:
                print(f'    - Descarga Img GOES *** UTC +00 *** {armado} - Ya existe en la carpeta.')
    
    return


  
def Procesar_GOES_MOD(selected_input_path, selected_output_path, utm = False, celsius = False, verbose = True):
    
    
    nc = rxr.open_rasterio(selected_input_path, mask_and_scale=True)

    # Se establece el la referencia espacial y se define en el archivo
    spatial_ref = 'PROJCS["unnamed",GEOGCS["unknown",DATUM["unnamed",SPHEROID["Spheroid",6378137,298.2572221]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]]],PROJECTION["Geostationary_Satellite"],PARAMETER["central_meridian",-75],PARAMETER["satellite_height",35786023],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],EXTENSION["PROJ4","+proj=geos +lon_0=-75 +h=35786023 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs +sweep=x"]]'
    nc_ref = nc.rio.set_crs(spatial_ref)

    # Transformación de coordenadas
    if utm == True:
        nc_ref = nc_ref.rio.reproject("EPSG:32720") 
    else:
        nc_ref = nc_ref.rio.reproject("EPSG:4326")

    # Clip por bbox
    nc_clip = nc_ref.rio.clip_box(
        minx=-68,
        miny=-42,
        maxx=-55,
        maxy=-21,
        crs="EPSG:4326",
    ).astype(float)

    nc_clip.rio.write_crs("EPSG:4326", inplace=True)

    # Cambio de Kelvin a Grados Celsius
    if celsius == True:
        nc_clip = nc_clip - 273.15

    # Se crea el nombre del archivo y dirección del archivo a guardar
    # file_name = sufijo_tif(file, source='goes')
    # tiff_name = f'{goes_wfolder}GOES_ABI-L2-LSTF_{file_name}'

    # Se guarda el .tif
    nc_clip.LST.rio.to_raster(selected_output_path)
    
    if(verbose):
            print(f'    - Preproceso GOES G1 - De .nc a .tiff. - Procesada Correctamente.')
            
    return

    

def Controlar_y_Procesar_GOES_G1(dir_root = "../", verbose = True):
    input_folder  = dir_root + "data/raw_data/original_data/GOES/"
    output_folder = dir_root + "data/raw_data/work_data/GOES/"

    input_files = os.listdir(input_folder)
    output_files = os.listdir(output_folder)

    pos_input = [23, 32]
    pos_output = [20, 29]

    jd_input = [x[pos_input[0]:pos_input[1]] for x in input_files]
    jd_output = [x[pos_output[0]:pos_output[1]] for x in output_files]
    
    # print(jd_input)
    # print(jd_output)
    dt_procesar = [x not in jd_output for x in jd_input]
    # print(dt_procesar)
    
    if sum(dt_procesar) > 0:
        
        if verbose:
            print(f'    - Preproceso GFS G1 - De .grib a .tif')
            print(f'    - Se tienen {sum(dt_procesar)} imagenes para preprocesar.')

        # Separamos los nombres de archivo de lo que hay que procesar
        # de input y el nombre que tendra en output
        input_file_a_procesar = [input_files[x] for x in range(len(input_files)) if dt_procesar[x]]
        output_file_nuevo = [NombrarGOES(x) for x in input_file_a_procesar]

        input_path_a_procesar  = [input_folder  + x for x in input_file_a_procesar] # path_input_a_procesar[x]
        output_path_a_procesar = [output_folder  + x for x in output_file_nuevo] # output_folder + output_file_nuevo[0]

        for x in range(len(input_path_a_procesar)):
            selected_input_path = input_path_a_procesar[x]
            selected_output_path = output_path_a_procesar[x]
            
            if verbose:
                print(f'    - Img {x + 1} de {len(input_path_a_procesar)}:')
                print(selected_input_path)
                print(selected_output_path)
                
            Procesar_GOES_MOD(selected_input_path, selected_output_path, utm = False, celsius = False, verbose = verbose)
            
            if verbose:
                print(" ")
                
            if x == range(len(input_path_a_procesar)):
                if verbose:
                    print(f'    - *** FIN DE Preproceso GFS G1 - De .grib a .tif ***')
    else:
        if verbose:
            print(f'    - Preproceso GFS G1 - De .grib a .tif. - Sin archivos para procesar - Todo OK!.')

    return


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # GFS
def nombre_gfs_tiff(file):

    # Toma el nombre del archivo .grib original del GFS
    fecha_jul = file.split("_")[-1][:-5] # corta el nombre por los '_' , se queda con el tercero contando desde el final. Corresponde a la fecha del inicio de la toma del dato.
    fecha_utc = pd.to_datetime(fecha_jul, format='%Y%m%d%H', utc=True) # se interpreta en formato de fecha en UTC

    # Guarda la fecha UTC en str en formato juliano y gregoriano
    j_utc = f'{fecha_utc:JAD%Y%j%H}' 
    g_utc = f'{fecha_utc:GAD%Y%m%d%H}' 

    # Convierte la fecha UTC a hora local
    fecha_local = fecha_utc.tz_convert('America/Cordoba')

    # Guarda la fecha local en str en formato juliano y gregoriano
    j_local = f'{fecha_local:JL%Y%j%H}'
    g_local = f'{fecha_local:GL%Y%m%d%H}'

    # Guarda en str el nombre completo del archivo tiff
    file_name = f'GFS_TSOIL_M00_{j_utc}_{g_utc}_{j_local}_{g_local}.tif'

    return file_name

############### Control para descarga GFS y GOES ###################
def lista_matcheo(path, source):
    ''' Genera un DataFrame con fecha y nombres
     de archivos .tif de una carpeta 
     
     Parámetro:
     path: ruta de la carpeta donde se encuentran los archivos
     '''
    if source == 'goes':
        ext = '.nc'
        ubi_date = -3
        prefix = 1
        suffix = 10
        format='%Y%j%H'
        path_files = [f for f in os.listdir(path) if f.endswith(ext)]
    elif source == 'gfs':
        ext = '.grib'
        ubi_date = -1
        prefix = 0
        suffix = 10
        format='%Y%m%d%H'
        path_files = [f for f in os.listdir(path) if f.endswith(ext)]
    elif source == 'prc':
        ext = '.tif'
        ubi_date = -3
        prefix = 3
        suffix = 13
        format='%Y%m%d%H'
        path_files = [f for f in os.listdir(path) if f.endswith(ext)]
    elif source == 'prc2':
        ext = '.tif'
        ubi_date = 3
        prefix = 3
        suffix = 13
        format='%Y%m%d%H'
        path_files = [f for f in os.listdir(path) if f.endswith(ext)]
        path_files = [x for x in path_files if 'GOES_ABI' in x]

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


def lista_procesar_proc(path_tif1, path_tif2, path_final):
    ''' 
    A partir de lista_matcheo genera dos listas con los nombres de archivos .tif cuyas fechas
    coniciden en dos carpetas distinas pero no en una tercera. Sirve para generar a partir de dos imágenes una tercera.

    Parámetros:
    path_tif1: Carpeta donde están los .tif a procesar
    path_tif2: Otra carpeta donde están los otros .tif a procesar
    path_final: Carpeta de destino de las imágenes procesadas
    '''

    path_tif1 = lista_matcheo(path_tif1, source='prc')
    path_tif2 = lista_matcheo(path_tif2, source='prc')
    path_final = lista_matcheo(path_final, source='prc2')

    procesar = pd.merge(path_tif1, path_tif2[["nom_archivo",'fecha_img']], left_on='fecha_img', right_on='fecha_img', how='left')
    procesar = pd.merge(procesar, path_final[["nom_archivo",'fecha_img']], left_on='fecha_img', right_on='fecha_img', how='left')
    procesar = procesar[pd.notnull(procesar.nom_archivo_y)]
    procesar = procesar[pd.isnull(procesar.nom_archivo)]
    
    list_procesar_x = list(procesar.nom_archivo_x)
    list_procesar_y = list(procesar.nom_archivo_y)

    return list_procesar_x, list_procesar_y


######## Control para descarga ##########
def control_descarga(start, end, path, source):
    ''' 
        Filtra una lista de fechas dada y entrega una lista de fechas que no han sido descargadas
        path: carpeta destino de las descargas
        source: fuente de las imágenes['goes' o 'gfs']
    '''
    if source == 'gfs':
        DATES = pd.date_range(start=start, end=end, freq="d")
        
        lista_horas = []
        for DATE in DATES:
            hora = 0
            for h_pred in range(24):
                d = f'{DATE:%Y%m%d}{hora:02d}'
                hora += 1
                lista_horas.append(d)

    elif source == 'goes':
        DATES = pd.date_range(start=start, end=end, freq="h")
        lista_horas = []
        for DATE in DATES:
            d = f'{DATE:%Y%m%d%H}'
            lista_horas.append(d)
    
    fechas_descarga = pd.DataFrame({'horas_descarga':lista_horas})
    
    en_carpeta = lista_matcheo(path,source=source)

    lista_descargar = pd.merge(fechas_descarga, en_carpeta[["fecha_img",'nom_archivo']], left_on='horas_descarga', right_on='fecha_img', how='left')
    lista_descargar = lista_descargar[pd.isnull(lista_descargar.nom_archivo)]
    lista_descargar = list(lista_descargar.horas_descarga)

    return lista_descargar



################# Descarga de GFS ######################
def descarga_gfs(start, end, dir_root = '../',verbose=True):
    warnings.filterwarnings("ignore")
    
    gfs_ofolder = f'{dir_root}data/raw_data/original_data/GFS/'
    
    lista_descarga = control_descarga(start=start,end=end,path=gfs_ofolder,source='gfs')

    if lista_descarga != []:

        # Se inerpretan las fechas
        start = pd.to_datetime(lista_descarga[0], format='%Y%m%d%H')
        end = pd.to_datetime(lista_descarga[-1], format='%Y%m%d%H')

        # Se genera una seria por hora entre las fechas seleccionadas
        DATES = pd.date_range(f"{start:%Y-%m-%d}", f"{end:%Y-%m-%d}", freq="D")    
        if verbose == True:
            print(f'** Inicio de la descarga GFS **')

        # Descarga el archivo que conincide con cada fecha
        for DATE in DATES:
            for h_pred in range(24):
                
                # Hora del modelo
                hora_mod = '00'

                # Se obtiene como str la hora predicha
                hora_pred = '{h:{fill}{width}}'.format(h=h_pred, fill='0', width=3)
                
                # URL de descarga de GFS TSOIL
                # url = (f'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t{hora_mod}z.pgrb2.0p25.f{hora_pred}&all_lev=on&var_TSOIL=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs.{DATE:%Y%m%d}%2F{hora_mod}%2Fatmos')
                
                # URL descarga de GFS TATM2m
                url = (f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t{hora_mod}z.pgrb2.0p25.f0{h_pred:02d}&lev_2_m_above_ground=on&var_TMP=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs.{DATE:%Y%m%d}%2F{hora_mod}%2Fatmos")

                # Se guarda el archivo original .grib
                file = requests.get(url)
                noaa_file = (f'gfs_0p25_m{hora_mod}_{DATE:%Y%m%d}{hora_pred[1:]}.grib')
                peso_archivo = open(gfs_ofolder+noaa_file, 'wb').write(file.content)
                
                # Si el archivo es chico, no concontró imagen y se borra el archivo generado
                if peso_archivo < 100000:
                    os.remove(gfs_ofolder+noaa_file)

                    if verbose == True:
                        print(f' [!] El modelo del {DATE:%Y-%m-%d} {hora_pred} UTC+0 no se encuentra disponible.')
                else:
                    
                    if verbose == True:
                        print(f'    - Imagen {DATE:%Y-%m-%d} {hora_pred} UTC+0 descargada correctamente.')
    else:
        if verbose == True:
            print('[!] No hay imagenes GFS nuevas para descargar')     
    return

# Función de procesamiento de GFS

def procesamiento_gfs(dir_root = '../', verbose=True):
    
    # Seteo de direcciones de apertura y guardado
    gfs_ofolder = f'{dir_root}data/raw_data/original_data/GFS/'
    gfs_wfolder = f'{dir_root}data/raw_data/work_data/GFS/'
    
    if verbose == True:
        print(f'** Inicio de procesamiento GFS **')
        if lista_procesar(gfs_ofolder, gfs_wfolder, source='gfs') == []:
            print('  - [!] No hay nuevas imagenes para procesar')    

    # Itera en todos los archivos encontrados en la carpeta original_data/WFS y procesa sólo que no tienen su correlativo procesado

    for file in lista_procesar(gfs_ofolder, gfs_wfolder, source='gfs'):
        
        pathfiles = f'{gfs_ofolder}{file}'

        gfs_grib = rxr.open_rasterio(pathfiles, mask_and_scale=True)
        gfs_band1 = xr.Dataset()
        gfs_band1['0-0.1[m] DBLL'] = gfs_grib[0]
        
        proj = CRS.from_wkt('GEOGCS["Coordinate System imported from GRIB file",DATUM["unnamed",SPHEROID["Sphere",6371229,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST]]')
        gfs_band1 = gfs_band1.rio.set_crs(proj)
        gfs_band1 = gfs_band1.rio.reproject("EPSG:4326")

        # Clip por bbox
        gfs_clip = gfs_band1.rio.clip_box(
            minx=-68,
            miny=-42,
            maxx=-55,
            maxy=-21,
            crs="EPSG:4326",
        )
        
        # Cambio de Grados Celsius a Kelvin
        gfs_clip = gfs_clip + 273.15
        
        # Se crea el nombre del archivo y dirección del archivo a guardar y se guarda
        file_name = sufijo_tif(file, source='gfs')
        tiff_name = f'{gfs_wfolder}GFS_TATM2_M00_{file_name}'
        gfs_clip.rio.to_raster(tiff_name)
        
        if verbose == True:
            ffile=file_name.split('_')[-3][3:]
            fecha_str = pd.to_datetime(ffile, format='%Y%m%d%H')
            print(f'    - Imagen {fecha_str:%Y-%m-%d %H:00} UTC+0 procesada correctamente.')

    return


def Detect_Descarga_Archivos_GFS(archivos_carpeta_gfs, archivos_online_gfs, dir_root = "../"):
    dt_descarga01 = [x in archivos_carpeta_gfs for x in archivos_online_gfs]
    dt_descarga01 = [x == False for x in dt_descarga01]
    return dt_descarga01


def Eliminar_1K_GFS_Files(one_gfs_path):
    
    eliminacion = False
    
    gfs_file_stats = os.stat(one_gfs_path)
    peso_archivo = gfs_file_stats.st_size
    limite = 1500
    if(peso_archivo < limite):
        os.remove(one_gfs_path) 
        eliminacion = True
    return eliminacion

def DescargarGfsMOD(dir_root = '../', verbose=True):
    
    # Harcodeado internamente para que base lo del dia 
    start = GFechaHora_Ingresa_Usuario_Descargar_GFS()
    
    if verbose == True:
        print(f'* Inicio de procesamiento GFS *')
        
    control_folder_gfs = f'{dir_root}data/raw_data/original_data/GFS/'
    archivos_carpeta_gfs = os.listdir(control_folder_gfs)
    
    # QUe ignore los avisos que surjan
    warnings.filterwarnings("ignore")

    # Guardamos la fecha gregoriana ingresada por el usuario
    # y la pasamos al formato intero
    gdf_original = start
    gdf_new = gdf_original[0:4] + gdf_original[5:7]  + gdf_original[8:10]
    # jdf_new = GLocal_to_JLocal(gdf_new)
    # jdf_original = jdf_new[x][0:4] + "/" + jdf_new[x][4:7] + "/" + jdf_new[x][7:10]
    
    # Hora del modelo
    hora_mod = '00'
    hora_total = 24
    horas_dia = range(hora_total)
    hora_pred = ['{h:{fill}{width}}'.format(h=x, fill='0', width=3) for x in horas_dia]

    general_url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t***hora_mod***z.pgrb2.0p25.f***hora_pred***&all_lev=on&var_TSOIL=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs.***gdf_new***%2F***hora_mod***%2Fatmos'
    general_file = 'gfs_0p25_m***hora_mod***_***gdf_new******hora_pred***.grib'
    all_url = []
    archivos_online_gfs = []
    for x in range(len(hora_pred)):
        nueva_url = general_url
        nueva_url = nueva_url.replace("***hora_mod***", hora_mod)
        nueva_url = nueva_url.replace("***hora_pred***", hora_pred[x])
        nueva_url = nueva_url.replace("***gdf_new***", gdf_new)
        all_url.append(nueva_url)
        nuevo_file = general_file
        nuevo_file = nuevo_file.replace("***hora_mod***", hora_mod)
        nuevo_file = nuevo_file.replace("***hora_pred***", hora_pred[x][1:])
        nuevo_file = nuevo_file.replace("***gdf_new***", gdf_new)
        archivos_online_gfs.append(nuevo_file)
        
        
    dt_descarga01 = Detect_Descarga_Archivos_GFS(archivos_carpeta_gfs, archivos_online_gfs, dir_root = "../")
    # dt_descarga01

    for x in range(len(dt_descarga01)):
        file_name = archivos_online_gfs[x]
        gd_file = file_name[13:23]
        jdf_file = GLocal_to_JLocal(gd_file)
        
        gd_user = gd_file[0:4] + "/" + gd_file[4:6] + "/" + gd_file[6:8] + " " + gd_file[8:10]
        jd_user = jdf_file[0:4] + "/" + jdf_file[4:7] + " " + jdf_file[7:9] 
        if dt_descarga01[x]:
            url_elegida = all_url[x]
            file_path = control_folder_gfs + file_name
            file = requests.get(url_elegida)


            open(file_path, 'wb').write(file.content)
            
            
            if verbose == True:
                print(f'    - Img GFS *** UTC +00 *** {jd_user} *** {gd_user} - Descargada Correctamente.')

        else:
            print(f'    - Img GFS *** UTC +00 *** {jd_user} *** {gd_user} - Ya esta en la carpeta.')



def Procesar_GFS_MOD(input_path, output_path, utm = False, celsius = False, verbose = True):
    
        
    gfs_grib = rxr.open_rasterio(input_path, mask_and_scale=True)
    gfs_band1 = xr.Dataset()
    gfs_band1['0-0.1[m] DBLL="Depth below land surface"'] = gfs_grib[0]

    # Detalle agregado por GERMAN.
    # De donde salio esto?
    proj = CRS.from_wkt('GEOGCS["Coordinate System imported from GRIB file",DATUM["unnamed",SPHEROID["Sphere",6371229,0]],PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST]]')

    #proj = CRS.from_epsg(9122) # se supone que esta imagen tiene sistema de coordenadas 9122
                               # como no viene con metadatos la imagen, se lo estamos imponiendo.

    gfs_band1 = gfs_band1.rio.set_crs(proj)
    gfs_band1 = gfs_band1.rio.reproject("EPSG:4326")

    # Clip por bbox
    gfs_clip = gfs_band1.rio.clip_box(
        minx=-68,
        miny=-42,
        maxx=-55,
        maxy=-21,
        crs="EPSG:4326",
    )

    # Cambio de Kelvin a Grados Celsius
    if celsius == True:
        gfs_clip = gfs_clip - 273.15 # cambio a celsius

    # Transformación de coordenadas
    if utm == True:
        gfs_clip = gfs_clip.rio.reproject("EPSG:32720")

    # Se crea el nombre del archivo y dirección del archivo a guardar y se guarda
    # file_name = nombre_gfs_tiff(file)
    # tiff_name = f'{gfs_wfolder}{file_name}'
    gfs_clip.astype(float).rio.to_raster(output_path)

    if(verbose):
        print(f'    - Preproceso GFS G1 - De .grib a .tiff. - Procesada Correctamente.')
    return


def Descarga_Controlada_GFS(dir_root = '../', start = GFechaHora_Ingresa_Usuario_Descargar_GFS(), verbose=True):
    
    # Harcodeado internamente para que base lo del dia    
    if verbose:
        print(f'* Inicio de descarga GFS *')
        
    control_folder_gfs = f'{dir_root}data/raw_data/original_data/GFS/'
    archivos_carpeta_gfs = os.listdir(control_folder_gfs)
    all_path_gfs = [control_folder_gfs + x for x in archivos_carpeta_gfs]
    
    dt_eliminados = []
    for este_archivo_gfs in all_path_gfs: 
        dt_eliminados.append(Eliminar_1K_GFS_Files(one_gfs_path = este_archivo_gfs))
    
    if sum(dt_eliminados) > 0: 
         if verbose:
                print(f'Se eliminaros {sum(dt_eliminados)} archivos GFS por ser archvios vacios.') 
    
    # Obtenemos nuevamente los archivos de la carpeta ya que pudimos haber eliminado archivos
    archivos_carpeta_gfs = os.listdir(control_folder_gfs)
    # QUe ignore los avisos que surjan
    warnings.filterwarnings("ignore")

    # Guardamos la fecha gregoriana ingresada por el usuario
    # y la pasamos al formato intero
    gdf_original = start
    gdf_new = gdf_original[0:4] + gdf_original[5:7]  + gdf_original[8:10]
    # jdf_new = GLocal_to_JLocal(gdf_new)
    # jdf_original = jdf_new[x][0:4] + "/" + jdf_new[x][4:7] + "/" + jdf_new[x][7:10]
    
    # Hora del modelo
    hora_mod = '00'
    hora_total = 24
    horas_dia = range(hora_total)
    hora_pred = ['{h:{fill}{width}}'.format(h=x, fill='0', width=3) for x in horas_dia]

    general_url = 'https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25_1hr.pl?file=gfs.t***hora_mod***z.pgrb2.0p25.f***hora_pred***&all_lev=on&var_TSOIL=on&leftlon=0&rightlon=360&toplat=90&bottomlat=-90&dir=%2Fgfs.***gdf_new***%2F***hora_mod***%2Fatmos'
    general_file = 'gfs_0p25_m***hora_mod***_***gdf_new******hora_pred***.grib'
    all_url = []
    archivos_online_gfs = []
    for x in range(len(hora_pred)):
        nueva_url = general_url
        nueva_url = nueva_url.replace("***hora_mod***", hora_mod)
        nueva_url = nueva_url.replace("***hora_pred***", hora_pred[x])
        nueva_url = nueva_url.replace("***gdf_new***", gdf_new)
        all_url.append(nueva_url)
        nuevo_file = general_file
        nuevo_file = nuevo_file.replace("***hora_mod***", hora_mod)
        nuevo_file = nuevo_file.replace("***hora_pred***", hora_pred[x][1:])
        nuevo_file = nuevo_file.replace("***gdf_new***", gdf_new)
        archivos_online_gfs.append(nuevo_file)
        
        
    dt_descarga01 = Detect_Descarga_Archivos_GFS(archivos_carpeta_gfs, archivos_online_gfs, dir_root = "../")
    # dt_descarga01

    for x in range(len(dt_descarga01)):
        file_name = archivos_online_gfs[x]
        gd_file = file_name[13:23]
        jdf_file = GLocal_to_JLocal(gd_file)
        
        gd_user = gd_file[0:4] + "/" + gd_file[4:6] + "/" + gd_file[6:8] + " " + gd_file[8:10]
        jd_user = jdf_file[0:4] + "/" + jdf_file[4:7] + " " + jdf_file[7:9] 
        if dt_descarga01[x]:
            url_elegida = all_url[x]
            file_path = control_folder_gfs + file_name
            file = requests.get(url_elegida)

            # Guardamos el archivo
            open(file_path, 'wb').write(file.content)
            
            # gfs_file_stats = os.stat(file)
            # peso_archivo = gfs_file_stats.st_size
            # print("El peso es ", peso_archivo)
            # Mandamos el sistema a dormir un toque
            time.sleep(1)
            
            dt_eliminacion = Eliminar_1K_GFS_Files(one_gfs_path = file_path)
            
            if verbose == True:
                if dt_eliminacion:
                    print(f'    - Img GFS *** UTC +00 *** {jd_user} *** {gd_user} - Archivo 1K - Descargado y Eliminado.')
                else:
                    print(f'    - Img GFS *** UTC +00 *** {jd_user} *** {gd_user} - Descargada Correctamente.')

        else:
            if verbose:
                print(f'    - Img GFS *** UTC +00 *** {jd_user} *** {gd_user} - Ya esta en la carpeta.')

    return

    
def Controlar_y_Procesar_GFS_G1(dir_root = "../", verbose = True):
    input_folder  = dir_root + "data/raw_data/original_data/GFS/"
    output_folder = dir_root + "data/raw_data/work_data/GFS/"

    input_files = os.listdir(input_folder)
    output_files = os.listdir(output_folder)

    pos_input = [13, 23]
    pos_output = [30, 40]

    jd_input = [x[pos_input[0]:pos_input[1]] for x in input_files]
    jd_output = [x[pos_output[0]:pos_output[1]] for x in output_files]

    dt_procesar = [x not in jd_output for x in jd_input]
    
    if sum(dt_procesar) > 0:
        if verbose:
            print(f'    - Preproceso GFS G1 - De .grib a .tif')
            print(f'    - Se tienen {sum(dt_procesar)} imagenes para preprocesar.')

        # Separamos los nombres de archivo de lo que hay que procesar
        # de input y el nombre que tendra en output
        input_file_a_procesar = [input_files[x] for x in range(len(input_files)) if dt_procesar[x]]
        # print(input_file_a_procesar)
        output_file_nuevo = [nombre_gfs_tiff(x) for x in input_file_a_procesar]

        input_path_a_procesar  = [input_folder  + x for x in input_file_a_procesar] # path_input_a_procesar[x]
        output_path_a_procesar = [output_folder  + x for x in output_file_nuevo] # output_folder + output_file_nuevo[0]
        # print(input_path_a_procesar)
        # print(output_path_a_procesar)
        
        for x in range(len(input_path_a_procesar)):
            selected_input_path = input_path_a_procesar[x]
            selected_output_path = output_path_a_procesar[x]
            
            if verbose:
                print(f'    - Img {x + 1} de {len(input_path_a_procesar)}:')
                print(selected_input_path)
                print(selected_output_path)
            
            Procesar_GFS_MOD(selected_input_path, selected_output_path, utm = False, celsius = False, verbose = verbose)
            
            if x == range(len(input_path_a_procesar)):
                if verbose:
                    print(f'    - *** FIN DE Preproceso GFS G1 - De .grib a .tif ***')
    else:
        if verbose:
            print(f'    - Preproceso GFS G1 - De .grib a .tif. - Sin archivos para procesar - Todo OK!.')

    return

##########################################################################################################################
##########################################################################################################################

def chequeo_internet():    
    # Se genera una seria por hora entre las fechas seleccionadas
    url = (f"http://www.google.com")
    try:
        response = requests.get(url, timeout=5)
        response.raise_for_status()
        internet = True
        # Code here will only run if the request is successful
    # except requests.exceptions.HTTPError as errh:
    #     internet = False
    except requests.exceptions.ConnectionError as errc:
        internet = False
    # except requests.exceptions.Timeout as errt:
    #     internet = False
    # except requests.exceptions.RequestException as err:
    #     internet = False
    return internet


##########################################################################################################################
##########################################################################################################################




##########################################################################################################################
##########################################################################################################################
