import sys
sys.path.append('../')
from funciones_G1 import *
from funciones_G2 import *
from funciones_G4 import *

import schedule
import time

def job():
    
    print("I'm working - Chequeo de conexión...")

    start = datetime.datetime.utcnow() - datetime.timedelta(days=1)
    end = datetime.datetime.utcnow()

    if chequeo_internet():
        
        verbose = True
        dir_root = '../'
        
        hora_inicio = datetime.datetime.now()
        print("Hora de Inicio Ejecucion: ", str(hora_inicio))
        
        print("I'm working - Generando carpetas...")
        crear_carpetas(dir_root=dir_root)
        
        print("I'm working - Descargando GOES de ayer...")
        Descarga_Controlada_GOES(dir_root = dir_root, gdh_usuario_descarga_goes = GFechaHora_descarga_utc_ayer_GOES(), verbose = verbose)
        # descarga_goes(dir_root=dir_root, start=start, end=end, verbose=verbose)

        print("I'm working - Descargando GOES de Hoy...")
        Descarga_Controlada_GOES(dir_root = dir_root, gdh_usuario_descarga_goes = GFechaHora_descarga_utc_hoy_GOES(), verbose = verbose)

        # print("I'm working - Descargando GFS de ayer...")
        # Descarga_Controlada_GFS(dir_root = dir_root, start = GFecha_ayerUTC(), verbose=verbose)
        # descarga_gfs(dir_root = dir_root, start = GFecha_ayerUTC(), verbose=verbose)

        print("I'm working - Descargando GFS de ayer y hoy...")
        # Descarga_Controlada_GFS(dir_root = dir_root, start = GFecha_hoyUTC(), verbose=verbose)
        descarga_gfs(dir_root=dir_root, start=start, end=end, verbose=verbose)

        print("I'm working - Procesando GOES...")
        Controlar_y_Procesar_GOES_G1(dir_root = dir_root, verbose = verbose)
        # procesamiento_goes(dir_root=dir_root, verbose=verbose)

        print("I'm working - Procesando GFS...")
        # Controlar_y_Procesar_GFS_G1(dir_root = "../", verbose = verbose)
        procesamiento_gfs(dir_root = dir_root, verbose = verbose)
        
        print("I'm working - Procesando Función Grupo 2...")
        proceso_grupo2(dir_root = dir_root, verbose = verbose)

        print("I'm working - Procesando Función Grupo 4 (Heladas acumuladas)...")
        hha(dir_root=dir_root, verbose = verbose)

        print(" ")
        hora_fin = datetime.datetime.now()
        print("Hora de Fin Ejecucion: ", str(hora_fin))

        tiempo_ejecucion = hora_fin - hora_inicio
        print("Tiempo de Ejecucion: ", str(tiempo_ejecucion))
        print("#############################################################")
    else:
        print('[!] No se puede acceder a la web para realizar las descargas')

schedule.every().hour.at(":12").do(job)

# python scripts/Cronos002.py

# .do(job)

# schedule.every().minute.at(":30").do(job)
# schedule.every().minute.at(":42").do(job)
# schedule.every(1).minutes.do(job)
# schedule.every(30).minutes.do(job)
# schedule.every(60).seconds.do(job)
while True:
    schedule.run_pending()
    time.sleep(1)