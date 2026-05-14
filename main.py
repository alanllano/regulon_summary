import os


def load_interactions(filename):
    """Lee interacciones del archivo TSV de reguladores y genes.

    Args:
        filename (str): Ruta al archivo TSV que contiene las interacciones.

    Returns:
        list[tuple[str, str, str]]: Lista de tuplas (TF, gene, effect), donde
            effect es uno de '+', '-', '-+'.
    """
    interactions = []
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
    return interactions


# "AraC" {
#    "genes": [araC, araA, araB, araD],
#    "activados": 4,
#    "reprimidos": 0
# }
def build_regulon(interactions):
    """Construye el regulon a partir de las interacciones leídas.

    Args:
        interactions (list[tuple[str, str, str]]): Lista de tuplas (TF, gene, effect).

    Returns:
        dict[str, dict[str, object]]: Diccionario con información por TF, incluyendo
            la lista de genes regulados y los conteos de activados y reprimidos.
    """
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
    return regulon


def get_tf_type(activados, reprimidos):
    """Determina el tipo de regulación del factor de transcripción.

    Args:
        activados (int): Número de genes activados por el TF.
        reprimidos (int): Número de genes reprimidos por el TF.

    Returns:
        str: 'Activador', 'Represor', 'Dual' o 'Desconocido'.
    """
    # Determinar del TF es activador, represor o dual
    if activados > 0 and reprimidos == 0:
        tipo_tf = "Activador"
    elif activados == 0 and reprimidos > 0:
        tipo_tf = "Represor"
    elif activados > 0 and reprimidos > 0:
        tipo_tf = "Dual"
    else:
        tipo_tf = "Desconocido"
    return tipo_tf


def write_summary(regulon, output_filename):
    """Escribe el resumen del regulon en un archivo TSV.

    Args:
        regulon (dict[str, dict[str, object]]): Información del regulon por TF.
        output_filename (str): Nombre del archivo de salida.
    """
    with open(output_filename, "w") as f:
        f.write(
            "TF\tTotal de genes regulados\tActivados\tReprimidos\tTipo TF\tLista de genes"
        )
        for tf in sorted(regulon):
            total_genes = len(regulon[tf]["genes"])
            lista_genes = ",".join(regulon[tf]["genes"])

            activados = regulon[tf]["activados"]
            reprimidos = regulon[tf]["reprimidos"]

            tipo_tf = get_tf_type(activados, reprimidos)

            f.write(
                f"{tf}\t{total_genes}\t{activados}\t{reprimidos}\t{tipo_tf}\t{lista_genes}\n"
            )


# =====
# main
# =====
def main():
    """Ejecuta el proceso completo de lectura, construcción y escritura del regulon.

    Lee el archivo de entrada, construye el regulon por TF y escribe un resumen
    en un archivo TSV de salida.
    """

    filename = "data/NetworkRegulatorGene.tsv"
    interactions = load_interactions(filename)

    regulon = build_regulon(interactions)

    output_filename = "regulon_summary.tsv"
    write_summary(regulon, output_filename)
    print(f"Archivo de salida generado: {output_filename}")


main()
