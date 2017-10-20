'''
Software for generating topographic surface and geoid model
with data from a global geoid model such as EGM08 and data 
for orthometic heights such as ETOPO2 model.
Copyright (C) 2017  Jose Ángel Hermosilla Rodrigo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from os import path, makedirs, SEEK_SET
from itertools import islice, chain
from functools import partial
from sys import argv, exit
from contextlib import ExitStack
from signal import signal, SIGINT
from shutil import rmtree
import math
import time

# Handler para el KeyboardInterrupt
def signal_handler(files, signal, frame):
    # Cerramos los archivos, de lo contrario no podremos eliminarlos
    for f in files: f.close()
    # Eliminamos los directorios enteros (con sus archivos)
    rmtree('salida_geoide_pm_{}'.format(paso_malla_salida))
    rmtree('salida_topo_pm_{}'.format(paso_malla_salida))
    # Salimos sin código de error
    exit(0)

# Códigos ANSI para los colores de la terminal
bcolors = {
    "HEADER"    : '\033[95m',
    "OKBLUE"    : '\033[94m',
    "OKGREEN"   : '\033[92m',
    "WARNING"   : '\033[93m',
    "FAIL"      : '\033[91m',
    "ENDC"      : '\033[0m' ,
    "BOLD"      : '\033[1m' ,
    "UNDERLINE" : '\033[4m' ,
}

# Función Helper para escribir en la terminal con los
# colores definidps en bcolors
def print_message(msg, color):
    print(bcolors[color] + msg + bcolors['ENDC'])

# Función para simplificar el archivo ETOPO
def simplifyETOPO(file):
    # Calculamos el paso de malla en un primer paso
    paso_malla_fichero = get_paso_malla_etopo(file)

    # Cada cuantas líneas debemos coger línea
    salto              = int(paso_malla_salida / paso_malla_fichero)

    # Índice que divide la lista en dos partes iguales
    middle             = math.floor( 360 * (60 / paso_malla_fichero) / 2 )

    def parse_line(l):

        # Pasamos la línea a floats
        line_values = list(map(float, l.strip().split()))

        # Como va de -180 a 180:
        # Generamos la lista de modo que queda [0...180, -180...0]
        # O lo que es lo mismo de [0...360]
        merged      = line_values[middle : -1] + line_values[ : middle]

        # Lo pasamos a metros y cogemos valores cada "salto"
        #return [ x*1000 for x in merged[::salto] ]
        return list( map( lambda x : x * 1000, merged[::salto] ) )
        
    #l = len([x for x in fileIn])
    # lines = []

    # for line in islice(file, 0, None, salto):
    #     lines += parse_line(line)
    # return lines

    # Un pelín más lento pero mola más
    return list( chain.from_iterable( map( parse_line, islice(file, 0, None, salto) ) ) )

# Función para simplificar el archivo EGM08
def simplifyEGM08(file):

    # Calculamos el paso de malla en un primer paso
    paso_malla_fichero = get_paso_malla_egm(file)

    # Este salto nos permite recoger las longitudes 
    # según el paso de malla de salida que queramos
    # Nos dice cuantas líneas hay en un grado de longitud
    salto              = int(paso_malla_salida / paso_malla_fichero)

    # Este salto nos permite saber hasta donde llega un
    # bloque de latitudes ej: Todas la lineas en que la latitud sea 89.5
    # Nos dice cuantas líneas hay para una misma latitud ej: Núm de línea para la lat 89.5
    salto_lat          = int(360 * 60 / paso_malla_fichero)

    # Número de líneas que hay que descartar tras haber leído un bloque
    # de latitudes. salto * salto_lat = Número de línes que hay en un grado de latitud
    # Le restamos un salto_lat (El que ya hemos leído)
    # Le sumamos 1 para coger la siguiente línea
    rest_of_lines      = salto * salto_lat - salto_lat + 1

    def parse_line(l):
        # Pasamos la lista a floats
        return list(map(float, l.strip().split()))

    # Lista que vamos a devolver
    lines = []

    # Bucle: Tenemos que repetir el bucle 181 veces
    # hasta llegar al final del archivo
    final = math.floor(180 * (60 / paso_malla_salida)) + 1
    for _ in range(0,  final):

        # Cogemos desde la linea "actual" hasta "actual + salto_lat_1", cada "salto"
        parsed_line = list(map(parse_line, islice(file, 0, salto_lat - 1, salto)))
        # Si es una lista vacía (Fin del archivo)

        #Añadimos la línea a la lista
        lines += parsed_line

        # Descartamos hasta el siguiente bloque
        for line in islice(file, rest_of_lines): pass

    return lines

# Convertir coordenadas geodésicas (lon, lat, h)
# en cartesianas geocéntricas (x, y, z)
def geo2car(lon, lat, h, factor_exageracion = 1):
    # Multiplicamos la h por la exageración
    h *= factor_exageracion
    
    # lon, lat en radianes
    lon_ = lon*math.pi/180
    lat_ = lat*math.pi/180

    # Semieje mayor WGS84
    a  = 6378137
    # Aplanamiento
    f  = 1.0/298.257223563
    # Semieje menor
    b  = a - (a*f)

    # excentricidad lineal
    el = math.sqrt((a**2) - (b**2))
    # primera excentricidad
    pe = el/a

    # v (nhu) radio primer vertical
    resta = ( (pe**2) * (math.sin(lat_)**2) )
    denom = math.sqrt( 1 - resta )
    v  = a / denom

    # Coordenadas finales x, y, z
    x  = (v + h)*math.cos(lat_)*math.cos(lon_)
    y  = (v + h)*math.cos(lat_)*math.sin(lon_)
    z  = ((b**2/a**2)*v + h)*math.sin(lat_)

    return [x, y, z]

# Convertir coordenadas cartesianas (x, y, z)
# a vértices conforme al estándar obj
def xyz2vertex(x, y, z):
    # Las dividimos entre 1.000.000 para 
    # obtener las coordenadas en unidades blender
    x = x / 1000000
    y = y / 1000000
    z = z / 1000000

    # Las formateamos con 6 decimales
    return 'v {:.6f} {:.6f} {:.6f}\n'.format(x, y, z)

# Obtener el paso de malla del fichero EGM
def get_paso_malla_egm(file):
    # Longitud inicial
    lon_i = None
    for line in file:
        # Obtenemos los valores de una línea
        lat_, lon_f, ond = list( map(float, line.strip().split()) )
        # Si la longitud inicial no está definida
        # La asignamos a la que hemos obtenido en esta primera línea
        if lon_i == None: lon_i = lon_f
        # Si no...
        else:
            # Obtenemos la diferencia de longitud entre la segunda y la primera línea
            dif = math.fabs(lon_f - lon_i)

            # Devolvemos el puntero de referencia del fichero a la posición 0 
            file.seek(SEEK_SET) # SEEK_SET es 0

            # Devolvemos el paso de malla calculado
            # Lo redondemos a dos decimales para evitarnos problemas
            # de decimales cuando hacemos operaciones de coma flotante
            return round(dif * 60, 2)

# Obtener el paso de malla del fichero ETOPO
def get_paso_malla_etopo(file):
    # Empezamos recorriendo el fichero
    # Con la primera línea nos basta
    # para saber el paso de malla
    for line in file:
        # Obtenemos los valores de la línea
        ls = line.strip().split()

        # Cuantos valores tiene la línea?
        num_elements = len(ls)

        # Devolvemos el puntero de referencia del fichero a la posición 0 
        file.seek(SEEK_SET) # SEEK_SET es 0

        # Devolvemos el paso de malla calculado
        return 60 / ((num_elements - 1 ) / 360)

if __name__ == '__main__':

    # Argumentos pasados por la línea de comandos
    args = argv[1:]
    # Tamaño de los argumentos
    args_len = len(args)

    # Si hay más de 4 ó 4 (Yoda Syntax)
    if 4 <= args_len :
        # Hemos pasado factores de exageración (Lista como números separados por coma)
        # Y el paso de malla de salida (int)
        filenameEGM08, filenameETOPO, factores_exageracion, paso_malla_salida = args[:4]
    # Hemos pasado 3 argumentos
    elif 3 == args_len:
        # Hemos pasado al menos los facotres de exageración
        filenameEGM08, filenameETOPO, factores_exageracion = args[:4]
        # Asignamos el paso de malla de salida por defecto
        print_message('Utilizando el paso de malla de salida por defecto: 60 minutos', 'WARNING')
        paso_malla_salida = 60
    # Si hemos pasado 2 argumentos
    elif 2 == args_len:
        # Hemos pasado al menos las rutas a los ficheros EGM y ETOPO
        filenameEGM08, filenameETOPO = args[:3]
        # Factores de exageración por defecto
        factores_exageracion = [1]
        print_message('Utilizando el factor de exageración por defecto: 1', 'WARNING')
        # Paso de malla de salida por defecto
        paso_malla_salida = 60
        print_message('Utilizando el paso de malla de salida por defecto: 60 minutos', 'WARNING')
    # Si no mostramos un error en rojo
    else:
        # Mostramos error
        print_message('No hay argumentos necesarios', 'FAIL')
        # Salimos del proceso
        exit(1)

    try:
        # Comprobamos que los ficheros existen
        # de lo contrario se lanzará un AssertionError
        # ya que estamos usando assert
        assert path.isfile(filenameEGM08)
        assert path.isfile(filenameETOPO)
        
        # Si hemos pasado factores de exageración, será un string
        # de lo contrario el valor por defecto que es [1]
        if not isinstance(factores_exageracion, list):
            # Generamos una lista con los valores pasados
            factores_exageracion = factores_exageracion.split(',')

        # Comprobamos que todos los elementos de la lista
        # pueden ser (o son) números enteros
        for x in factores_exageracion:
            # No es entero, será str
            if not type(x) == int:
                # No es un dígito (i.e: Valor no convertible a int)
                if not x.isdigit():
                    # Lanzamos excepción
                    raise ValueError('Factor de exageración debe contener número enteros')

        # Convertimos a lista de enteros
        # para cercionarnos de que tenemos int y no string en la lista
        factores_exageracion = list(map(int, factores_exageracion))

        # print(factores_exageracion)

        # Si el paso de malla de salida es de tipo str
        # es el que habremos pasado por la línea de comandos
        if type(paso_malla_salida) == str:
            # Si no representa un número entero
            if not paso_malla_salida.isdigit():
                # Lanzamos excepción
                raise ValueError('Paso de malla debe ser un número entero')
            # Convertimos a int el paso de malla de salida
            paso_malla_salida  = int(paso_malla_salida)
    # Manejamos excepciones
    except AssertionError as e:
        print_message('Ruta/s no válida/s', 'FAIL')
        exit(1)
    except ValueError as e:
        print_message(e.__str__(), 'FAIL')
        exit(1)

    # Nombre de los directorios de salida que vamos a crear
    # Uno para el geoide, otro para la sup. topográfica
    outdir_geoide = 'salida_geoide_pm_{}'.format(paso_malla_salida)
    outdir_topo = 'salida_topo_pm_{}'.format(paso_malla_salida)

    # Creamos los directorios de salida si no existen
    makedirs(outdir_geoide, exist_ok = True)
    makedirs(outdir_topo, exist_ok = True)

    print_message('Carpetas de salida creadas', 'OKGREEN')

    # Abrimos los ficheros fileEGM y fileTopo
    # y un objeto de la clase ExitStack
    # que crea un Stack que nos permite abrir
    # ficheros, de forma programática y cerrarlos automáticamente
    # al terminar el estamento with
    with ExitStack() as stack, \
        open(filenameEGM08, 'rt') as fileEGM, \
        open(filenameETOPO, 'rt') as fileTopo:

        # Vamos a abrir todos los ficheros de escritura a la vez
        # Uno por cada factor de exageración
        filenames_geoide = map(lambda fe: path.join(outdir_geoide, 
                                    'geoide_f{}_pm{}.obj'.format(fe, paso_malla_salida)), \
                        factores_exageracion)

        filenames_topo = map(lambda fe: path.join(outdir_topo, 
                                    'topo_f{}_pm{}.obj'.format(fe, paso_malla_salida)), \
                        factores_exageracion)

        # Abrimos o Creamos los archivos de escritura
        files_geoide = [ stack.enter_context(open(fname, 'wt')) for fname in filenames_geoide ]
        files_topo = [ stack.enter_context(open(fname, 'wt')) for fname in filenames_topo ]

        # Señal que se dispara cuando ocurre un KeyBoardInterrupt
        files = list(files_geoide) + list(files_topo)
        signal(SIGINT, partial(signal_handler, files))

        # Tiempo
        t = time.process_time()

        egm  = simplifyEGM08(fileEGM)
        topo = simplifyETOPO(fileTopo)

        # Loop sobre el fihcero egm
        for i, l in enumerate(egm):
            # Valores del fichero EGM
            lat, lon, ond = egm[i]
            # Valor del fichero etopo
            hort          = topo[i]

            # Calculamos la altura elipsoidal
            helip         = hort + ond

            # Para cada factor de exageración ...
            for idx, fe in enumerate(factores_exageracion):
                # Obtenemos los ficheros en los que debemos escribir
                # según el factor de exageración
                file_geoide = files_geoide[idx]
                file_topo = files_topo[idx]

                # Escribimos los resultados de obtener las cordenadas de la sup. topográfica
                # y el geoide respectivamente para el factor de exageración actual
                file_geoide.write( xyz2vertex( *geo2car(lon, lat, ond, fe) ) )
                file_topo.write( xyz2vertex( *geo2car( lon, lat, helip, fe) ) )

        # Tiempo que ha tardado
        elapsed = time.process_time() - t
        print_message('Ha tardado {}s'.format(elapsed), 'OKGREEN')