import os

# =========================================
# Lectura del archivo y construcción de interactions
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

            print(f"Leído: {line[:50]}...")

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
            if len(fields) <= 7:
                continue

            # columnas a utilizar
            TF = fields[1]
            gene = fields[3]
            effect = fields[4]

            # Validar effect
            if effect not in ["+", "-"]:
                continue

            interactions.append((TF, gene, effect))

print(interactions)
