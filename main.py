# Leer el archivo de entrada
# Contar # genes activados y reprimidos para cada TF
# Generar archivo de salida
# Calcular si un TF es activador o represor o dual

# Refactorizar el código para que sea más legible y eficiente


import os

# =========================================
# Lectura del archivo y construcción de interactions
# Responsabilidad: leer el archivo de entrada y construir una lista de tuplas con la información relevante
# Entrada: archivo tsv con información de reguladores y genes
# Salida: lista de tuplas (TF, gene, effect)
# =========================================

interactions = []

filename = "data/raw/NetworkRegulatorGene.tsv"

if not os.path.exists(filename):
    print("Error: archivo no encontrado")
    exit(1)
else:
    with open(filename) as f:

        for line in f:

            line = line.strip()

            # Ignorar líneas vacías
            if not line:
                continue

            # Ignorar comentarios
            if line.startswith("#"):
                continue

            # Ignorar encabezado
            if line.startswith("1)regulatorId"):
                continue

            fields = line.split("\t")

            # Validar número mínimo de columnas
            if len(fields) <= 6:
                continue

            # columnas a utilizar
            TF = fields[1]
            gene = fields[4]
            effect = fields[5]

            # Validar effect
            if effect not in ["+", "-", "-+"]:
                continue

            interactions.append((TF, gene, effect))


# =========================================
# Construcción del regulon con información extra
# Entrada: lista de tuplas (TF, gene, effect)
# Salida: diccionario con TF como clave y lista de genes, conteo de activados y reprimidos como valor
# =========================================
regulon = {}  # diccionario con lista de genes
for tf, gene, effect in interactions:
    if tf not in regulon:
        regulon[tf] = {"genes": [], "activados": 0, "reprimidos": 0}
    regulon[tf]["genes"].append(gene)

    # Contar activados y reprimidos
    if effect == "+":
        regulon[tf]["activados"] += 1
    elif effect == "-":
        regulon[tf]["reprimidos"] += 1
    elif effect == "-+":
        regulon[tf]["activados"] += 1
        regulon[tf]["reprimidos"] += 1

# "AraC" {
#    "genes": [araC, araA, araB, araD],
#    "activados": 4,
#    "reprimidos": 0
# }


# =========================================
# Generación de la salida
# Entrada: diccionario con TF como clave y lista de genes, conteo de activados y reprimidos como valor
# Salida: archivo tsv con resumen de cada TF
# imprimir en un archivo el resumen de cada TF
# =========================================
output_filename = "regulon_summary.tsv"
with open(output_filename, "w") as f:
    f.write(
        "TF\tTotal de genes regulados\tActivados\tReprimidos\tTipo TF\tLista de genes"
    )
    for tf in sorted(regulon):
        total_genes = len(regulon[tf]["genes"])
        lista_genes = ",".join(regulon[tf]["genes"])
        activados = regulon[tf]["activados"]
        reprimidos = regulon[tf]["reprimidos"]

        # Determinar del TF es activador, represor o dual
        if activados > 0 and reprimidos == 0:
            tipo_tf = "Activador"
        elif activados == 0 and reprimidos > 0:
            tipo_tf = "Represor"
        elif activados > 0 and reprimidos > 0:
            tipo_tf = "Dual"

        f.write(
            f"{tf}\t{total_genes}\t{activados}\t{reprimidos}\t{tipo_tf}\t{lista_genes}\n"
        )
print(f"Archivo de salida generado: {output_filename}")
