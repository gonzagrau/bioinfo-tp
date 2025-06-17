import os
from pathlib import Path

def combinar_fastas():
    # Definir rutas
    carpeta_fastas = "outputs/emboss/orf_splits"
    archivo_salida = "orf_splits_combined.fasta"
    
    # Verificar si la carpeta existe
    if not os.path.exists(carpeta_fastas):
        print(f"Error: No se encontró la carpeta '{carpeta_fastas}'")
        return
    
    # Obtener lista de archivos FASTA en la carpeta
    archivos_fasta = [f for f in os.listdir(carpeta_fastas) if f.endswith(('.fasta', '.fa', '.fna', '.ffn'))]

    if not archivos_fasta:
        print(f"No se encontraron archivos FASTA en la carpeta '{carpeta_fastas}'")
        return
    
    # Crear archivo combinado
    with open(archivo_salida, 'w') as outfile:
        for archivo in archivos_fasta:
            ruta_archivo = os.path.join(carpeta_fastas, archivo)
            print(f"Procesando: {archivo}")
            
            
            # Copiar el contenido del archivo FASTA
            with open(ruta_archivo, 'r') as infile:
                for linea in infile:
                    outfile.write(linea)
    
    print(f"\n¡Combinación completada! Resultado guardado en '{archivo_salida}'")
    print(f"Se combinaron {len(archivos_fasta)} archivos FASTA.")

if __name__ == "__main__":
    combinar_fastas()