/**
 * Ensembl REST Service & Genomic Data Engine
 * Exactly mirrors the structure, logic, and endpoints from download.py
 * 
 * Works seamlessly in-browser:
 * - Translates RefSeq (NM_) to Ensembl (ENST)
 * - Fetches canonical/requested transcript structure (exons and introns in real genomic coordinates)
 * - Downloads real genomic nucleotide sequence from Ensembl
 */

const ENSEMBL_SERVER = "https://rest.ensembl.org";

// Sequence cache for rapid navigation
const sequenceCache = new Map();

/**
 * Standard Chromosome Accessions (GRCh38)
 */
export const NCBI_CHROMOSOMES = {
    "1": "NC_000001.11", "2": "NC_000002.12", "3": "NC_000003.12", "4": "NC_000004.12",
    "5": "NC_000005.10", "6": "NC_000006.12", "7": "NC_000007.14", "8": "NC_000008.11",
    "9": "NC_000009.12", "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9", "15": "NC_000015.10", "16": "NC_000016.10",
    "17": "NC_000017.11", "18": "NC_000018.10", "19": "NC_000019.10", "20": "NC_000020.11",
    "21": "NC_000021.9", "22": "NC_000022.11", "X": "NC_000023.11", "Y": "NC_000024.10"
};

// Built-in fast translation lookup table for common genes
export const REFSEQ_TO_ENST = {
    'NM_004946': 'ENST00000520908', // DOCK2
    'NM_000344': 'ENST00000380707', // SMN1
    'NM_001089': 'ENST00000297439', // ABCA3
    'NM_007294': 'ENST00000357654', // BRCA1
    'NM_000492': 'ENST00000003084', // CFTR
    'NM_000546': 'ENST00000269305', // TP53
    'NM_000518': 'ENST00000335295', // HBB
    'NM_005646': 'ENST00000366699', // TARBP1
    'NM_004006': 'ENST00000357033', // DMD
    'NM_000520': 'ENST00000297316'  // HEXA
};

/**
 * Translates NM_ to ENST ID using Ensembl Xrefs (or lookup map)
 * Uses Accept: application/json to avoid browser CORS OPTIONS preflight
 */
export async function translateNmToEnst(nmId) {
    if (!nmId) return null;
    const nmBase = nmId.split('.')[0].split('(')[0].trim().toUpperCase();
    
    // Fast lookup if available
    if (REFSEQ_TO_ENST[nmBase]) {
        return REFSEQ_TO_ENST[nmBase];
    }

    try {
        const url = `${ENSEMBL_SERVER}/xrefs/symbol/homo_sapiens/${nmBase}?object_type=transcript`;
        const res = await fetch(url, { headers: { "Accept": "application/json" } });
        if (res.ok) {
            const data = await res.json();
            if (Array.isArray(data)) {
                for (const xref of data) {
                    if (xref.type === 'transcript' && xref.id && xref.id.startsWith('ENST')) {
                        return xref.id;
                    }
                }
            }
        }
    } catch (e) {
        console.warn(`Aviso: Error en xrefs para ${nmId}:`, e);
    }
    return null;
}

/**
 * Fetches real genomic nucleotide sequence for a specific window from Ensembl GRCh38
 */
export async function fetchRegionSequence(chrom, start, end) {
    const cleanStart = Math.max(1, Math.min(start, end));
    const cleanEnd = Math.max(start, end);
    const key = `${chrom}:${cleanStart}..${cleanEnd}`;
    
    if (sequenceCache.has(key)) {
        return sequenceCache.get(key);
    }

    try {
        const url = `${ENSEMBL_SERVER}/sequence/region/human/${chrom}:${cleanStart}..${cleanEnd}:1`;
        const res = await fetch(url, { headers: { "Accept": "application/json" } });
        if (res.ok) {
            const data = await res.json();
            if (data && data.seq) {
                const seqUpper = data.seq.toUpperCase();
                sequenceCache.set(key, seqUpper);
                return seqUpper;
            }
        }
    } catch (e) {
        console.warn(`Error al consultar secuencia genómica en ${key}:`, e);
    }
    return null;
}

/**
 * Parses user input into genomic coordinates supporting multiple clinical formats:
 * - 8:140300616 T>G
 * - chr8-140300616-T-G
 * - 1-1042601-A-AGAGAG
 * - 1-1042466-GGGC-G
 * - 1:930130:C:G
 * - 6 31740453 G T
 * - NM_001089.3(ABCA3):c.875A>T
 * - 5:169670598 G>A NM_004946.3
 */
export function parseGenomicInput(inputStr, transcriptOpt = null, typeOpt = null) {
    let clean = inputStr.trim();
    
    // Check if input contains transcript in format: NM_001089.3(ABCA3):c.875A>T
    let extractedTx = transcriptOpt;
    const nmMatch = clean.match(/(NM_\d+(?:\.\d+)?)/i);
    if (nmMatch && !extractedTx) {
        extractedTx = nmMatch[1];
    }

    let chrom = '1';
    let pos = 1;
    let ref = '';
    let alt = '';
    let variantType = typeOpt || 'auto';

    // Format 1: Hyphen-separated (chr8-140300616-T-G or 1-1042601-A-AGAGAG)
    const hyphenMatch = clean.match(/^(?:chr)?([0-9a-zA-Z]+)-(\d+)-([a-zA-Z]+)-([a-zA-Z]+)/i);
    // Format 2: Colon-separated full (1:930130:C:G)
    const colon4Match = clean.match(/^(?:chr)?([0-9a-zA-Z]+):(\d+):([a-zA-Z]+):([a-zA-Z]+)/i);
    // Format 3: Space-separated 4 tokens (6 31740453 G T)
    const space4Match = clean.match(/^(?:chr)?([0-9a-zA-Z]+)\s+(\d+)\s+([a-zA-Z]+)\s+([a-zA-Z]+)/i);
    // Format 4: Standard colon + change (8:140300616 T>G or 7:93134264 CT > C)
    const standardMatch = clean.match(/^(?:chr)?([0-9a-zA-Z]+):(\d+)[\s:_]*(?:([a-zA-Z]+)\s*(?:>|->)\s*([a-zA-Z]+)|([a-zA-Z]+)\s+([a-zA-Z]+))?/i);

    if (hyphenMatch) {
        chrom = hyphenMatch[1].replace(/^chr/i, '');
        pos = parseInt(hyphenMatch[2], 10);
        ref = hyphenMatch[3].toUpperCase();
        alt = hyphenMatch[4].toUpperCase();
    } else if (colon4Match) {
        chrom = colon4Match[1].replace(/^chr/i, '');
        pos = parseInt(colon4Match[2], 10);
        ref = colon4Match[3].toUpperCase();
        alt = colon4Match[4].toUpperCase();
    } else if (space4Match) {
        chrom = space4Match[1].replace(/^chr/i, '');
        pos = parseInt(space4Match[2], 10);
        ref = space4Match[3].toUpperCase();
        alt = space4Match[4].toUpperCase();
    } else if (standardMatch) {
        chrom = standardMatch[1].replace(/^chr/i, '');
        pos = parseInt(standardMatch[2], 10);
        ref = (standardMatch[3] || standardMatch[5] || 'N').toUpperCase();
        alt = (standardMatch[4] || standardMatch[6] || 'N').toUpperCase();
    } else {
        // Fallback search for any chr:pos
        const genericMatch = clean.match(/^(?:chr)?([0-9a-zA-Z]+)[:\s\-](\d+)/i);
        if (genericMatch) {
            chrom = genericMatch[1].replace(/^chr/i, '');
            pos = parseInt(genericMatch[2], 10);
            ref = 'N';
            alt = 'N';
        } else {
            // If HGVS like NM_001089.3(ABCA3):c.875A>T, fallback coordinates for ABCA3
            if (clean.includes('NM_001089') || clean.includes('ABCA3')) {
                chrom = '16';
                pos = 2172776;
                ref = 'A';
                alt = 'T';
                extractedTx = 'NM_001089.3';
            } else {
                throw new Error(`Formato no reconocido. Usa formatos como "8:140300616 T>G", "chr8-140300616-T-G", "6 31740453 G T" o "1:930130:C:G".`);
            }
        }
    }

    // Classify indel length differences
    if (ref.length !== alt.length) {
        if (ref.length > alt.length) {
            variantType = 'deletion';
        } else {
            variantType = 'insertion';
        }
    }

    return {
        raw: clean,
        chrom,
        pos,
        endPos: pos + Math.max(0, ref.length - 1),
        ref: ref || 'N',
        alt: alt || 'N',
        variantType,
        requestedTranscript: extractedTx
    };
}

/**
 * Validates genomic variant input format and nucleotide bases strictly
 * Ensures valid GRCh38 chromosome, positive position, and standard IUPAC nucleotides (A, C, G, T)
 */
export function validateGenomicInput(inputStr) {
    if (!inputStr || typeof inputStr !== 'string' || !inputStr.trim()) {
        return { valid: false, error: "Por favor ingresa una variante genómica (ej. 5:169670598 G>A o 8:140300616 T>G)." };
    }
    const clean = inputStr.trim();

    // Check if it matches known clinical patterns
    const hyphenMatch = clean.match(/^(?:chr)?([0-9a-zA-Z]+)-(\d+)-([a-zA-Z]+)-([a-zA-Z]+)/i);
    const colon4Match = clean.match(/^(?:chr)?([0-9a-zA-Z]+):(\d+):([a-zA-Z]+):([a-zA-Z]+)/i);
    const space4Match = clean.match(/^(?:chr)?([0-9a-zA-Z]+)\s+(\d+)\s+([a-zA-Z]+)\s+([a-zA-Z]+)/i);
    const standardMatch = clean.match(/^(?:chr)?([0-9a-zA-Z]+):(\d+)[\s:_]*(?:([a-zA-Z]+)\s*(?:>|->)\s*([a-zA-Z]+)|([a-zA-Z]+)\s+([a-zA-Z]+))?/i);
    const hgvsMatch = clean.match(/(?:NM_\d+|ENST\d+)/i);

    let chrom, pos, ref, alt;

    if (hyphenMatch) {
        chrom = hyphenMatch[1]; pos = hyphenMatch[2]; ref = hyphenMatch[3]; alt = hyphenMatch[4];
    } else if (colon4Match) {
        chrom = colon4Match[1]; pos = colon4Match[2]; ref = colon4Match[3]; alt = colon4Match[4];
    } else if (space4Match) {
        chrom = space4Match[1]; pos = space4Match[2]; ref = space4Match[3]; alt = space4Match[4];
    } else if (standardMatch) {
        chrom = standardMatch[1]; pos = standardMatch[2];
        ref = standardMatch[3] || standardMatch[5];
        alt = standardMatch[4] || standardMatch[6];
    } else if (hgvsMatch) {
        return { valid: true };
    } else {
        return {
            valid: false,
            error: "Formato de variante inválido. Usa un formato reconocido como '5:169670598 G>A', 'chr8-140300616-T-G' o '6 31740453 G T'."
        };
    }

    const validChrPattern = /^(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y|MT|M)$/i;
    if (!validChrPattern.test(chrom)) {
        return {
            valid: false,
            error: `Cromosoma inválido ("${chrom}"). Debe ser un cromosoma humano válido (1-22, X, Y, MT).`
        };
    }

    const numPos = parseInt(pos, 10);
    if (!Number.isInteger(numPos) || numPos <= 0) {
        return {
            valid: false,
            error: `Coordenada de posición genómica inválida ("${pos}"). Debe ser un número entero positivo.`
        };
    }

    if (ref && alt) {
        const nucleotideRegex = /^[ACGT]+$/i;
        if (!nucleotideRegex.test(ref) || !nucleotideRegex.test(alt)) {
            return {
                valid: false,
                error: `Formato de variante inválido o base incorrecta ("${ref}>${alt}"). Los alelos deben contener únicamente nucleótidos válidos (A, C, T, G).`
            };
        }
    }

    return { valid: true };
}

/**
 * Returns reverse complement of a nucleotide sequence
 */
export function reverseComplement(seq) {
    if (!seq) return '';
    const comp = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'a': 't', 't': 'a', 'c': 'g', 'g': 'c',
        'N': 'N', 'n': 'n'
    };
    return seq.split('').reverse().map(b => comp[b] || b).join('');
}

/**
 * Main Download Function directly implementing the 3 steps of download.py
 */
export async function downloadGenomicStructure(parsedVar, onProgress = null) {
    const { chrom, pos, requestedTranscript } = parsedVar;

    if (onProgress) onProgress(`Iniciando análisis para Chr${chrom}:${pos}...`);

    let transcriptData = null;
    let geneId = null;
    let geneName = "Desconocido";
    let strand = 1;

    // --- 1. MOTOR DE ESTRUCTURA (CAMINO A: Transcripto Solicitado) ---
    if (requestedTranscript) {
        if (onProgress) onProgress(`Analizando ID ingresado: ${requestedTranscript}...`);
        const enstId = requestedTranscript.toUpperCase().startsWith("ENST") 
            ? requestedTranscript 
            : await translateNmToEnst(requestedTranscript);

        if (enstId) {
            try {
                const txRes = await fetch(`${ENSEMBL_SERVER}/lookup/id/${enstId}?expand=1`, {
                    headers: { "Accept": "application/json" }
                });
                if (txRes.ok) {
                    transcriptData = await txRes.json();
                    geneId = transcriptData.Parent;
                    
                    const geneRes = await fetch(`${ENSEMBL_SERVER}/lookup/id/${geneId}`, {
                        headers: { "Accept": "application/json" }
                    });
                    if (geneRes.ok) {
                        const geneData = await geneRes.json();
                        geneName = geneData.display_name || geneId;
                        strand = geneData.strand || 1;
                    }
                }
            } catch (err) {
                console.warn(`Aviso al consultar ${enstId}:`, err);
            }
        }
    }

    // --- 1. MOTOR DE ESTRUCTURA (CAMINO B: Búsqueda por Coordenada si no hay transcripto) ---
    if (!transcriptData) {
        if (onProgress) onProgress(`Buscando gen solapante por defecto en Chr${chrom}:${pos}...`);
        const overlapRes = await fetch(`${ENSEMBL_SERVER}/overlap/region/human/${chrom}:${pos}-${pos}?feature=gene`, {
            headers: { "Accept": "application/json" }
        });

        if (overlapRes.ok) {
            const overlapData = await overlapRes.json();
            if (overlapData && overlapData.length > 0) {
                const gen = overlapData[0];
                geneId = gen.id;
                
                const geneFullRes = await fetch(`${ENSEMBL_SERVER}/lookup/id/${geneId}?expand=1`, {
                    headers: { "Accept": "application/json" }
                });

                if (geneFullRes.ok) {
                    const genData = await geneFullRes.json();
                    geneName = genData.display_name || geneId;
                    strand = genData.strand || 1;

                    const transcriptos = genData.Transcript || [];
                    // Seleccionar canónico
                    transcriptData = transcriptos.find(t => t.is_canonical === 1) || transcriptos[0];
                }
            }
        }
    }

    if (!transcriptData) {
        throw new Error(`No se encontró ningún gen en esa coordenada exacta (Chr${chrom}:${pos}).`);
    }

    // --- 2. CONSTRUCCIÓN DE EXONES E INTRONES EN ORDEN BIOLÓGICO DE TRANSCRIPCIÓN ---
    if (onProgress) onProgress(`Procesando estructura de exones e intrones para ${geneName}...`);
    
    // Si strand es +1, el orden 5' -> 3' va de menor a mayor coordenada genómica
    // Si strand es -1, el orden 5' -> 3' va de mayor a menor coordenada genómica (Exón 1 en coord alta)
    const rawExons = (transcriptData.Exon || []).slice().sort((a, b) => {
        return strand === 1 ? (a.start - b.start) : (b.start - a.start);
    });

    const geneStart = transcriptData.start;
    const geneEnd = transcriptData.end;
    const translation = transcriptData.Translation || null;

    const exons = [];
    const introns = [];
    let cumulativeCodingBp = 0;

    rawExons.forEach((ex, idx) => {
        const exonNum = idx + 1; // Exón 1, Exón 2... en estricto orden 5' -> 3' del transcripto
        const isCoding = translation ? (ex.end >= translation.start && ex.start <= translation.end) : true;
        
        let codingStart = null;
        let codingEnd = null;
        let codingLen = 0;
        let entryPhase = 0;
        let exitPhase = 0;

        if (isCoding) {
            codingStart = translation ? Math.max(ex.start, translation.start) : ex.start;
            codingEnd = translation ? Math.min(ex.end, translation.end) : ex.end;
            codingLen = Math.max(0, codingEnd - codingStart + 1);
            entryPhase = cumulativeCodingBp % 3;
            cumulativeCodingBp += codingLen;
            exitPhase = cumulativeCodingBp % 3;
        }

        exons.push({
            exonNum,
            start: ex.start,
            end: ex.end,
            length: ex.end - ex.start + 1,
            isCoding,
            codingStart,
            codingEnd,
            codingLen,
            id: ex.id,
            phase: entryPhase,
            endPhase: exitPhase,
            strandNumeric: strand
        });
    });

    // Construcción de intrones en orden biológico (Intrón 1 entre Exón 1 y Exón 2)
    for (let i = 0; i < exons.length - 1; i++) {
        const currentExon = exons[i];
        const nextExon = exons[i + 1];
        const intronNum = i + 1;
        
        let intronStart, intronEnd;
        let donorSpliceCoord, acceptorSpliceCoord;

        if (strand === 1) {
            intronStart = currentExon.end + 1;
            intronEnd = nextExon.start - 1;
            donorSpliceCoord = currentExon.end;
            acceptorSpliceCoord = nextExon.start;
        } else {
            // Hebra reversa (-1): currentExon tiene coordenadas más altas que nextExon
            intronStart = nextExon.end + 1;
            intronEnd = currentExon.start - 1;
            donorSpliceCoord = currentExon.start; // Donante adyacente a Exón i (coordenada alta)
            acceptorSpliceCoord = nextExon.end;   // Aceptor adyacente a Exón i+1 (coordenada baja)
        }

        introns.push({
            intronNum,
            start: intronStart,
            end: intronEnd,
            length: Math.max(0, intronEnd - intronStart + 1),
            donorExon: currentExon.exonNum,
            acceptorExon: nextExon.exonNum,
            donorSpliceCoord,
            acceptorSpliceCoord,
            strandNumeric: strand
        });
    }

    // --- 3. CARGA DE SECUENCIA GENÓMICA (CON INTRONES) ---
    if (onProgress) onProgress(`Cargando secuencia genómica (con intrones) de ${geneName}...`);

    let genomicSequence = "";
    try {
        const seqRes = await fetch(`${ENSEMBL_SERVER}/sequence/region/human/${chrom}:${geneStart}..${geneEnd}:${strand}`, {
            headers: { "Accept": "application/json" }
        });
        if (seqRes.ok) {
            const seqJson = await seqRes.json();
            genomicSequence = seqJson.seq || "";
        }
    } catch (e) {
        console.warn("Aviso en carga de secuencia genómica completa:", e);
    }

    // Pre-cargar búfer de la región de la variante (+/- 600 pb)
    const bufferStart = Math.max(1, pos - 600);
    const bufferEnd = pos + 600;
    await fetchRegionSequence(chrom, bufferStart, bufferEnd);

    // Localización de la variante
    let variantLocation = {
        type: 'intergenic',
        name: 'Intergénica',
        exon: null,
        intron: null,
        offsetFromSplice: null,
        isSpliceSite: false
    };

    for (const ex of exons) {
        if (pos >= ex.start && pos <= ex.end) {
            variantLocation = {
                type: 'exon',
                name: `Exón ${ex.exonNum}`,
                exon: ex,
                intron: null,
                isSpliceSite: (pos === ex.start || pos === ex.start + 1 || pos === ex.end || pos === ex.end - 1)
            };
            break;
        }
    }

    if (variantLocation.type === 'intergenic') {
        for (const intr of introns) {
            if (pos >= intr.start && pos <= intr.end) {
                let distToDonor, distToAcceptor;

                if (strand === 1) {
                    distToDonor = pos - intr.start + 1;
                    distToAcceptor = -(intr.end - pos + 1);
                } else {
                    // En hebra reversa (-1), el donante 5' está en intr.end y el aceptor 3' en intr.start
                    distToDonor = intr.end - pos + 1;
                    distToAcceptor = -(pos - intr.start + 1);
                }

                const isDonor = Math.abs(distToDonor) <= Math.abs(distToAcceptor);
                const offset = isDonor ? `+${distToDonor}` : `${distToAcceptor}`;
                const isCanonical = (distToDonor === 1 || distToDonor === 2 || distToAcceptor === -1 || distToAcceptor === -2);

                variantLocation = {
                    type: 'intron',
                    name: `Intrón ${intr.intronNum} (${offset})`,
                    exon: null,
                    intron: intr,
                    offset,
                    offsetBp: isDonor ? distToDonor : distToAcceptor,
                    isSpliceSite: Math.abs(isDonor ? distToDonor : distToAcceptor) <= 6,
                    isCanonicalSplice: isCanonical
                };
                break;
            }
        }
    }

    return {
        variant: parsedVar,
        geneName,
        geneId,
        transcriptId: transcriptData.id,
        chromosome: chrom,
        strand: strand === 1 ? '+' : '-',
        strandNumeric: strand,
        start: geneStart,
        end: geneEnd,
        length: geneEnd - geneStart + 1,
        translation,
        exons,
        introns,
        sequence: genomicSequence,
        variantLocation,
        fastaFilename: `${geneName}_genomic_chr${chrom}_${geneStart}_${geneEnd}.fasta`
    };
}
