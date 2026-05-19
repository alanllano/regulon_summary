import argparse
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

    with open(filename, encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            # Ignorar líneas vacías
            if not line:
                continue

            # Ignorar comentarios
            if line.startswith("#"):
                continue

            # Ignorar encabezado de columnas si existe
            if line.startswith("1)riId") or line.startswith("1)regulatorId"):
                continue

            fields = line.split("\t")

            # Validar número mínimo de columnas
            if len(fields) < 19:
                continue

            # columnas a utilizar según el formato de TF-RISet
            TF = fields[3].strip()
            gene = fields[16].strip() or fields[18].strip()
            effect_raw = fields[10].strip().lower()

            if not TF or not gene or not effect_raw:
                continue

            if effect_raw in ["+", "activator", "activation", "positive", "up"]:
                effect = "+"
            elif effect_raw in ["-", "repressor", "repression", "negative", "down"]:
                effect = "-"
            elif effect_raw in [
                "-+",
                "dual",
                "dual regulator",
                "mixed",
                "dual-function",
            ]:
                effect = "-+"
            else:
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
        if gene not in regulon[tf]["genes"]:
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
    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(
            "TF\tTotal de genes regulados\tActivados\tReprimidos\tTipo TF\tLista de genes\n"
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
# argumentos
# =====
def read_arguments():
    parser = argparse.ArgumentParser(
        description="Construye un regulon a partir de un archivo TSV de interacciones entre TF y genes."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Archivo TSV de entrada con interacciones entre TF y genes",
    )
    parser.add_argument("-o", "--output", required=True)
    return parser.parse_args()


# =====
# main
# =====
def main():
    """Ejecuta el proceso completo de lectura, construcción y escritura del regulon.

    Lee el archivo de entrada, construye el regulon por TF y escribe un resumen
    en un archivo TSV de salida.
    """

    args = read_arguments()

    filename = args.input
    interactions = load_interactions(filename)

    regulon = build_regulon(interactions)

    output_filename = args.output

    write_summary(regulon, output_filename)
    print(f"Archivo de salida generado: {output_filename}")


if __name__ == "__main__":
    main()
