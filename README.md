# 🧬 Aiuda-predictor — Genomic VariantStudio

**Aiuda-predictor** es una plataforma bioinformática interactiva y pedagógica para la consulta de estructura génica (exones e intrones en coordenadas GRCh38), visualización a nivel macro y micro de secuencias nucleotídicas, y simulación en tiempo real del impacto de variantes de splicing y exónicas (desplazamiento de marco de lectura / frameshift, codones de parada prematuros PTC y activación de degradación mediada por NMD).

---

## 🚀 Características Principales

1. **Paso 1: Descarga y Mapeo Genómico Automático (Ensembl REST API):**
   - Resolución de variantes en múltiples formatos clínicos (`Chr:Pos Ref>Alt`, `Chr-Pos-Ref-Alt`, `NM_001089.3(ABCA3):c.875A>T`, etc.).
   - Mapeo automático de transcritos canónicos y descarga de estructura exón-intrón real en GRCh38.
   - Cálculo estricto de fases de lectura de entrada y salida (`Fase 0 ➔ 2`, `Fase 1 ➔ 0`, etc.) para cada exón codificante.

2. **Paso 2: Mapa Macro de Exones e Intrones (Visualiza tu variante en el gen):**
   - Vista balanceada y vista a escala genómica real.
   - Centrado automático interactivo en el locus de la variante.
   - Tabla de exones plegada por default con navegación rápida a los extremos de cada elemento.

3. **Paso 3: Visor a Nivel de Bases y Límite de Splicing (A nivel de bases):**
   - Visualizador de nucleótidos (coloreados por base) con regla genómica, demarcación de límites de empalme y visualización de codones/aminoácidos en chevrons.
   - Soporte completo para variantes puntuales (SNVs), inserciones y deleciones.
   - Conmutador entre **Secuencia Salvaje (WT)** y **Variante Mutada**.

4. **Paso 4: Clasificación y Simulación de la Variante (Predecí qué pasará):**
   - **Clasificación del Tipo de Variante:** Deducción del impacto biológico y clínico.
   - **Simulador de Splicing (SpliceAI + 3 Escenarios):** Exclusivo para variantes de splicing canónico o cercanas a las uniones exón-intrón.
     - *Salto de Exón Completo (Exon Skipping):* Comparación 1:1 en 3 filas escalonadas.
     - *Activación de Sitios Crípticos (Cryptic Sites):* Estructura primaria con corte y el exón aceptor a la derecha, seguido del nuevo empalme con el aminoácido híbrido de unión.
     - *Retención de Intrón Completo (Intron Retention):* Búsqueda en secuencia real hasta el primer codón STOP prematuro.
   - **Simulador de Marco de Lectura Exónico:** Recálculo de marco de lectura para variantes exónicas.

---

## 💻 Instalación y Uso Local

No requiere configuración compleja ni librerías externas en el servidor. Funciona nativamente con tecnologías web modernas (HTML5, Vanilla CSS y ES Modules).

### 1. Clonar el repositorio:
```bash
git clone https://github.com/EGCinmuno/Aiuda-predictor.git
cd Aiuda-predictor
```

### 2. Levantar el servidor local:
Puedes usar Python:
```bash
python -m http.server 8080
```
O Node.js con cualquier servidor estático:
```bash
npx serve .
```

### 3. Abrir en el navegador:
Accede a `http://localhost:8080`.

---

## 🛠️ Tecnologías

- **Frontend:** HTML5 Semántico, Vanilla CSS (Variables CSS, Diseño responsivo, Modo oscuro/claro, Glassmorphism).
- **Lógica & Bioinformática:** JavaScript ES Modules (Vanilla JS nativo).
- **APIs:** [Ensembl REST API](https://rest.ensembl.org/) (Genoma de Referencia Humano GRCh38).

---

## 📄 Licencia

Desarrollado para **EGC**.
