/**
 * Visual Renderer Engine
 * 
 * 1. renderGenomicExonMap() — Non-overlapping, balanced & linear exon-intron map with Inicio/Final jump buttons & loading states
 * 2. renderFranklinBaseViewer() — Micro base-level viewer with dynamic region headers [ Exón 4 ], prominent focus banner, and fixed +/-10pb navigation
 * 3. renderSynchronizedGridTrack() — Synchronized 1:1 column base & amino acid track (with pre-mRNA untranslated mode & strict STOP truncation)
 * 4. renderComparisonSplicingViewer() — Full pedagogical index-based WT vs Splicing viewers:
 *    - Caso 1: Exon Skipping ([Final Exón anterior] -> Se perdería: [Inicio saltado]--//--[Final saltado] -> [Inicio siguiente] -> [Resultado: Empalme final])
 *    - Caso 2: Cryptic Site (1. Estructura Primaria Pre-ARNm, 2. Empalme Mutado con regiones identificadas, mutación visible y traducción)
 *    - Caso 3: Retención de Intrón (Marco continuo hasta el STOP real en secuencia genómica, corte estricto)
 * 5. renderConsequenceSimulator() — Step 6 interactive consequence evaluator
 */

import { CODON_TABLE, aa1to3, aaFullName } from './geneticCode.js';
import { fetchRegionSequence, reverseComplement } from './ensemblService.js';
import { fetchSpliceAIPrediction } from './spliceaiService.js';

// Base colors
export const BASE_COLORS = {
    'A': '#059669', // Emerald Green
    'C': '#0284c7', // Sky Blue
    'G': '#f59e0b', // Amber/Orange
    'T': '#f43f5e', // Rose Red
    'N': '#64748b'  // Slate
};

// Phase colors for exons
export const PHASE_COLORS = {
    0: '#0284c7',   // Cyan/Blue (In-frame 3n)
    1: '#8b5cf6',   // Purple (Split codon 1 nt)
    2: '#10b981',   // Emerald (Split codon 2 nt)
    null: '#475569' // Slate UTR
};

/**
 * Cache for split codon lookups across splice junctions
 */
const splitCodonCache = new Map();

export async function getSplitCodonInfo(currentExon, exons, chromosome) {
    if (!currentExon || !currentExon.isCoding || currentExon.endPhase === 0) {
        return null;
    }
    const isAntisense = (currentExon.strandNumeric === -1);
    const cacheKey = `${chromosome}_E${currentExon.exonNum}_end_${isAntisense ? 'neg' : 'pos'}`;
    if (splitCodonCache.has(cacheKey)) {
        return splitCodonCache.get(cacheKey);
    }

    const splitCount1 = currentExon.endPhase; // 1 or 2
    const splitCount2 = 3 - splitCount1;      // 2 or 1
    
    const nextExon = exons.find(e => e.exonNum === currentExon.exonNum + 1 && e.isCoding);
    if (!nextExon) return null;

    try {
        let tailSeq, headSeq;
        if (!isAntisense) {
            tailSeq = await fetchRegionSequence(chromosome, currentExon.end - splitCount1 + 1, currentExon.end);
            headSeq = await fetchRegionSequence(chromosome, nextExon.start, nextExon.start + splitCount2 - 1);
        } else {
            // Hebra reversa (-1): 3' end del currentExon en currentExon.start, 5' de nextExon en nextExon.end
            const rawTail = await fetchRegionSequence(chromosome, currentExon.start, currentExon.start + splitCount1 - 1);
            const rawHead = await fetchRegionSequence(chromosome, nextExon.end - splitCount2 + 1, nextExon.end);
            tailSeq = reverseComplement(rawTail);
            headSeq = reverseComplement(rawHead);
        }

        if (tailSeq && headSeq) {
            const fullCodon = (tailSeq + headSeq).toUpperCase();
            const aa = CODON_TABLE[fullCodon] || 'X';
            const aa3 = aa1to3(aa);
            const name = aaFullName(aa);

            const info = {
                splitCount1,
                splitCount2,
                tailSeq,
                headSeq,
                fullCodon,
                aa,
                aa3,
                name,
                currentExonNum: currentExon.exonNum,
                nextExonNum: nextExon.exonNum
            };
            splitCodonCache.set(cacheKey, info);
            return info;
        }
    } catch (e) {
        console.warn("Aviso en split codon lookup:", e);
    }
    return null;
}

/**
 * Helper to get a codon-aligned tail segment of an exon in 5' -> 3' biological orientation.
 * Ensures sequence slices in Step 6 start on clean codon boundaries (Phase 0).
 */
export function getCodonAlignedExonTail(exon, targetLen = 15) {
    if (!exon || !exon.isCoding) {
        return { start: exon ? exon.start : 1, end: exon ? exon.end : 1, phase: 0 };
    }
    const isAntisense = (exon.strandNumeric === -1);
    if (isAntisense) {
        // En hebra (-1), el 3' tail biológico está en exon.start; el 5' en exon.end
        let p = Math.min(exon.end, exon.start + targetLen - 1);
        const rem = (exon.end - p + exon.phase) % 3;
        if (rem !== 0) {
            if (p + rem <= exon.end) p += rem;
            else p -= (3 - rem);
        }
        return { start: exon.start, end: p, phase: 0 };
    } else {
        let pos = Math.max(exon.start, exon.end - targetLen + 1);
        const rem = (pos - exon.start + exon.phase) % 3;
        if (rem !== 0) {
            if (pos - rem >= exon.start) {
                pos = pos - rem; // retroceder a límite de codón
            } else {
                pos = pos + (3 - rem); // avanzar a límite de codón
            }
        }
        return {
            start: pos,
            end: exon.end,
            phase: 0 // Fase 0 garantizada
        };
    }
}

/**
 * Helper to get an exon head segment in 5' -> 3' biological orientation starting with native phase.
 */
export function getCodonAlignedExonHead(exon, targetLen = 15) {
    if (!exon) return { start: 1, end: 1, phase: 0 };
    const isAntisense = (exon.strandNumeric === -1);
    if (isAntisense) {
        // En hebra (-1), el 5' head biológico está en exon.end
        return {
            start: Math.max(exon.start, exon.end - targetLen + 1),
            end: exon.end,
            phase: exon.phase ?? 0
        };
    } else {
        return {
            start: exon.start,
            end: Math.min(exon.end, exon.start + targetLen - 1),
            phase: exon.phase ?? 0
        };
    }
}

/**
 * Fetches an exon segment and formats it in 5' -> 3' mRNA orientation
 */
export async function fetchExonSegment5to3(chromosome, segInfo, strandNumeric) {
    let raw = await fetchRegionSequence(chromosome, segInfo.start, segInfo.end);
    if (!raw) raw = "N".repeat(segInfo.end - segInfo.start + 1);
    return strandNumeric === -1 ? reverseComplement(raw) : raw;
}

/**
 * ═══════════════════════════════════════════════════════════════
 * PASO 3: MAPA GENERAL DE EXONES E INTRONES (SIN SOLAPAMIENTO)
 * ═══════════════════════════════════════════════════════════════
 */
let exonMapViewMode = 'balanced'; // 'balanced' or 'linear'

export function renderGenomicExonMap(container, model, onRegionClick = null) {
    if (!container || !model) return;
    container.innerHTML = '';

    const { exons, introns, start, end, strand, chromosome, variant, variantLocation } = model;
    const totalGenomicSpan = Math.max(1, end - start);

    const wrapper = document.createElement('div');
    wrapper.className = 'exon-map-wrapper';

    // Top Controls Bar
    const controlsBar = document.createElement('div');
    controlsBar.className = 'exon-map-controls-bar';
    controlsBar.innerHTML = `
        <div class="exon-count-badge">
            📊 <strong>${exons.length} Exones</strong> &bull; <strong>${introns.length} Intrones</strong> &bull; Hebra (${strand})
        </div>
        <div style="display: flex; gap: 8px; align-items: center; flex-wrap: wrap;">
            <div class="map-view-toggle-group">
                <button id="btnViewBalanced" class="btn-map-view ${exonMapViewMode === 'balanced' ? 'active' : ''}">
                    ✨ Vista Estructura Balanceada
                </button>
                <button id="btnViewLinear" class="btn-map-view ${exonMapViewMode === 'linear' ? 'active' : ''}">
                    📏 Vista Escala Genómica Real
                </button>
            </div>
            <button id="toggleExonTableBtn" class="btn-secondary" style="padding: 6px 14px; font-size: 0.8rem; cursor: pointer;">
                🔽 Plegar Tabla de Exones
            </button>
        </div>
    `;
    wrapper.appendChild(controlsBar);

    const trackWrap = document.createElement('div');
    trackWrap.className = 'exon-map-track-container';

    if (exonMapViewMode === 'balanced') {
        const diagramTrack = document.createElement('div');
        diagramTrack.className = 'exon-diagram-track';

        exons.forEach((ex) => {
            const isVariantExon = (variantLocation.type === 'exon' && variantLocation.exon?.exonNum === ex.exonNum);
            const intronAfter = introns.find(intr => intr.donorExon === ex.exonNum || intr.intronNum === ex.exonNum);
            const isVariantIntron = (variantLocation.type === 'intron' && variantLocation.intron?.intronNum === intronAfter?.intronNum);

            const exonCard = document.createElement('div');
            exonCard.className = `exon-diagram-card ${isVariantExon ? 'is-variant-target' : ''}`;
            exonCard.style.borderColor = isVariantExon ? 'var(--accent-rose)' : (PHASE_COLORS[ex.phase] || '#0284c7');

            let phaseBadge = ex.isCoding 
                ? `<span class="exon-card-phase phase-${ex.phase}">Fase ${ex.phase} ➔ ${ex.endPhase}</span>`
                : `<span class="exon-card-phase">UTR</span>`;

            exonCard.innerHTML = `
                ${isVariantExon ? `
                    <div class="diagram-variant-pin">
                        <span>📍 ${variant.ref}&gt;${variant.alt}</span>
                    </div>
                ` : ''}
                <div class="exon-card-header">
                    <span class="exon-card-num">Exón ${ex.exonNum}</span>
                    ${phaseBadge}
                </div>
                <div class="exon-card-length">${ex.length.toLocaleString()} pb</div>
                <div class="exon-card-coords">Chr${chromosome}:${ex.start.toLocaleString()}–${ex.end.toLocaleString()}</div>
            `;

            exonCard.addEventListener('click', () => {
                if (onRegionClick) onRegionClick(ex.start);
            });
            diagramTrack.appendChild(exonCard);

            if (intronAfter) {
                const intronWrap = document.createElement('div');
                intronWrap.className = `intron-connector-wrap ${isVariantIntron ? 'is-variant-target' : ''}`;
                
                const intronKb = (intronAfter.length / 1000).toFixed(1);
                intronWrap.innerHTML = `
                    ${isVariantIntron ? `
                        <div class="diagram-variant-pin intron-pin">
                            <span>📍 ${variantLocation.offset} (${variant.ref}&gt;${variant.alt})</span>
                        </div>
                    ` : ''}
                    <div class="intron-connector-line">
                        <span class="intron-flow-arrow">${strand === '+' ? '►' : '◄'}</span>
                    </div>
                    <div class="intron-connector-label" title="Intrón ${intronAfter.intronNum} (${intronAfter.length.toLocaleString()} pb)">
                        Intrón ${intronAfter.intronNum} (${intronKb >= 1 ? `${intronKb} kb` : `${intronAfter.length} pb`})
                    </div>
                `;

                intronWrap.addEventListener('click', () => {
                    if (onRegionClick) onRegionClick(intronAfter.start);
                });
                diagramTrack.appendChild(intronWrap);
            }
        });

        trackWrap.appendChild(diagramTrack);
    } else {
        const linearTrack = document.createElement('div');
        linearTrack.className = 'exon-map-track';

        const backbone = document.createElement('div');
        backbone.className = 'exon-backbone';
        linearTrack.appendChild(backbone);

        const strandMarker = document.createElement('div');
        strandMarker.className = 'strand-label';
        strandMarker.textContent = strand === '+' ? '5′ ──────► 3′' : '3′ ◄────── 5′';
        linearTrack.appendChild(strandMarker);

        exons.forEach(ex => {
            const leftPct = ((ex.start - start) / totalGenomicSpan) * 100;
            const widthPct = Math.max(0.6, ((ex.end - ex.start + 1) / totalGenomicSpan) * 100);

            const exonEl = document.createElement('div');
            exonEl.className = 'exon-block';
            exonEl.style.left = `${leftPct}%`;
            exonEl.style.width = `${widthPct}%`;

            const isVariantInThisExon = (variantLocation.type === 'exon' && variantLocation.exon?.exonNum === ex.exonNum);
            exonEl.style.background = isVariantInThisExon 
                ? 'var(--accent-rose)' 
                : (ex.isCoding ? (PHASE_COLORS[ex.phase] || PHASE_COLORS[0]) : PHASE_COLORS[null]);

            exonEl.title = `Exón ${ex.exonNum} | Chr${chromosome}:${ex.start.toLocaleString()} - ${ex.end.toLocaleString()} (${ex.length} pb)`;

            if (widthPct > 2.0) {
                const lbl = document.createElement('span');
                lbl.className = 'exon-label';
                lbl.textContent = `E${ex.exonNum}`;
                exonEl.appendChild(lbl);
            }

            exonEl.addEventListener('click', () => {
                if (onRegionClick) onRegionClick(ex.start);
            });

            linearTrack.appendChild(exonEl);
        });

        if (variant && variant.pos) {
            const varPct = Math.max(0, Math.min(100, ((variant.pos - start) / totalGenomicSpan) * 100));
            const pin = document.createElement('div');
            pin.className = 'variant-lollipop';
            pin.style.left = `${varPct}%`;
            pin.innerHTML = `
                <div class="lollipop-head" title="Posición genómica: Chr${chromosome}:${variant.pos.toLocaleString()}">
                    📍 ${variant.ref} &gt; ${variant.alt}
                </div>
                <div class="lollipop-stem"></div>
            `;
            linearTrack.appendChild(pin);
        }

        trackWrap.appendChild(linearTrack);
    }

    wrapper.appendChild(trackWrap);

    // Detailed Table
    const tableWrap = document.createElement('div');
    tableWrap.id = 'exonDetailsTableContainer';
    tableWrap.className = 'exon-table';
    tableWrap.innerHTML = `
        <table>
            <thead>
                <tr>
                    <th>Elemento</th>
                    <th>Inicio Genómico</th>
                    <th>Fin Genómico</th>
                    <th>Tamaño</th>
                    <th>Fase de Lectura</th>
                    <th>Navegación al Visor (Paso 4)</th>
                </tr>
            </thead>
            <tbody>
                ${exons.map((ex) => {
                    const isExonHit = (variantLocation.type === 'exon' && variantLocation.exon?.exonNum === ex.exonNum);
                    const intronAfter = introns.find(intr => intr.donorExon === ex.exonNum || intr.intronNum === ex.exonNum);
                    const isIntronHit = (variantLocation.type === 'intron' && variantLocation.intron?.intronNum === intronAfter?.intronNum);

                    let phaseBadgeHtml = '<span class="phase-badge">UTR</span>';
                    if (ex.isCoding) {
                        const entryDesc = ex.phase === 0 ? 'Codón completo al inicio' : ex.phase === 1 ? '1 nt previo (split 2/3 al inicio)' : '2 nt previos (split 1/3 al inicio)';
                        const exitDesc = ex.endPhase === 0 ? 'Codón termina exacto al final' : ex.endPhase === 1 ? '1 nt remanente [1/3 ➔]' : '2 nt remanentes [2/3 ➔]';
                        phaseBadgeHtml = `<span class="phase-badge phase-${ex.phase}" title="Inicio: Fase ${ex.phase} (${entryDesc}). Final: Fase ${ex.endPhase} (${exitDesc}).">Fase ${ex.phase} ➔ ${ex.endPhase}</span>`;
                    }

                    let rows = `
                        <tr class="${isExonHit ? 'variant-exon-row' : ''}">
                            <td><strong>Exón ${ex.exonNum}</strong></td>
                            <td>${ex.start.toLocaleString()}</td>
                            <td>${ex.end.toLocaleString()}</td>
                            <td>${ex.length.toLocaleString()} pb</td>
                            <td>${phaseBadgeHtml}</td>
                            <td>
                                <div style="display: flex; gap: 4px; align-items: center; flex-wrap: wrap;">
                                    <button class="btn-secondary jump-exon-btn" data-coord="${ex.start}" title="Centrar visor en el inicio de Exón ${ex.exonNum} (Chr${chromosome}:${ex.start})" style="padding: 3px 8px; font-size: 0.72rem; font-weight: 600;">
                                        🟢 Ver Inicio
                                    </button>
                                    <button class="btn-secondary jump-exon-btn" data-coord="${ex.end}" title="Centrar visor en el final de Exón ${ex.exonNum} (Chr${chromosome}:${ex.end})" style="padding: 3px 8px; font-size: 0.72rem; font-weight: 600;">
                                        🔴 Ver Final
                                    </button>
                                    ${isExonHit ? `<span class="badge badge-danger" style="margin-left: 4px;">📍 Variante</span>` : ''}
                                </div>
                            </td>
                        </tr>
                    `;

                    if (intronAfter) {
                        rows += `
                            <tr class="${isIntronHit ? 'variant-exon-row' : ''}" style="opacity: 0.85;">
                                <td style="padding-left: 20px; color: var(--text-muted);">↳ <em>Intrón ${intronAfter.intronNum}</em></td>
                                <td>${intronAfter.start.toLocaleString()}</td>
                                <td>${intronAfter.end.toLocaleString()}</td>
                                <td>${intronAfter.length.toLocaleString()} pb</td>
                                <td><span class="phase-badge" style="background: rgba(6,182,212,0.15); color: #38bdf8; border: 1px solid rgba(6,182,212,0.4);">Intrón</span></td>
                                <td>
                                    <div style="display: flex; gap: 4px; align-items: center; flex-wrap: wrap;">
                                        <button class="btn-secondary jump-exon-btn" data-coord="${intronAfter.start}" title="Inicio de Intrón ${intronAfter.intronNum} (sitio dador +1)" style="padding: 3px 8px; font-size: 0.72rem;">
                                            🔍 Inicio (+1)
                                        </button>
                                        <button class="btn-secondary jump-exon-btn" data-coord="${intronAfter.end}" title="Final de Intrón ${intronAfter.intronNum} (sitio aceptor -1)" style="padding: 3px 8px; font-size: 0.72rem;">
                                            🔍 Final (-1)
                                        </button>
                                        ${isIntronHit ? `<span class="badge badge-warning" style="margin-left: 4px;">📍 Splicing ${variantLocation.offset}</span>` : ''}
                                    </div>
                                </td>
                            </tr>
                        `;
                    }
                    return rows;
                }).join('')}
            </tbody>
        </table>
    `;
    wrapper.appendChild(tableWrap);

    container.appendChild(wrapper);

    // Auto-center on variant element in Step 3
    setTimeout(() => {
        const targetElement = wrapper.querySelector('.exon-diagram-card.is-variant-target, .intron-connector-wrap.is-variant-target, .variant-lollipop');
        if (targetElement) {
            targetElement.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'center' });
        }
    }, 60);

    // Event Listeners for view mode toggle
    const btnBalanced = wrapper.querySelector('#btnViewBalanced');
    const btnLinear = wrapper.querySelector('#btnViewLinear');
    if (btnBalanced && btnLinear) {
        btnBalanced.addEventListener('click', () => {
            exonMapViewMode = 'balanced';
            renderGenomicExonMap(container, model, onRegionClick);
        });
        btnLinear.addEventListener('click', () => {
            exonMapViewMode = 'linear';
            renderGenomicExonMap(container, model, onRegionClick);
        });
    }

    // Toggle table collapse/expand
    const toggleBtn = wrapper.querySelector('#toggleExonTableBtn');
    if (toggleBtn && tableWrap) {
        toggleBtn.addEventListener('click', () => {
            const isCollapsed = tableWrap.style.display === 'none';
            tableWrap.style.display = isCollapsed ? 'block' : 'none';
            toggleBtn.textContent = isCollapsed ? '🔽 Plegar Tabla de Exones' : '▶ Desplegar Tabla de Exones';
        });
    }

    // Jump to exon buttons with Loading state feedback
    wrapper.querySelectorAll('.jump-exon-btn').forEach(btn => {
        btn.addEventListener('click', async () => {
            const coord = parseInt(btn.dataset.coord, 10);
            if (onRegionClick && !isNaN(coord)) {
                const origHtml = btn.innerHTML;
                btn.disabled = true;
                btn.innerHTML = '<span class="calc-spinner">⏳ Cargando secuencia...</span>';
                try {
                    await onRegionClick(coord);
                } finally {
                    btn.disabled = false;
                    btn.innerHTML = origHtml;
                }
            }
        });
    });
}

/**
 * ═══════════════════════════════════════════════════════════════
 * PASO 4: VISOR A NIVEL DE BASES Y UNIÓN DE SPLICING
 * ═══════════════════════════════════════════════════════════════
 */
const BASE_WINDOW_SIZE = 30;
let isShowingMutantView = false;

export async function renderFranklinBaseViewer(container, model, customCenterPos = null, showMutant = null) {
    if (!container || !model) return;

    if (showMutant !== null) {
        isShowingMutantView = showMutant;
    }

    try {
        const { exons, introns, chromosome, variant, variantLocation } = model;
        const isAntisense = (model.strandNumeric === -1 || model.strand === '-');
        
        let targetPos = variant.pos;
        if (customCenterPos !== null && customCenterPos !== undefined) {
            targetPos = customCenterPos;
        }

        const half = Math.floor(BASE_WINDOW_SIZE / 2);
        let winStart = Math.max(1, targetPos - half);
        let winEnd = winStart + BASE_WINDOW_SIZE - 1;
        const windowLength = winEnd - winStart + 1;

        let realSequence = await fetchRegionSequence(chromosome, winStart, winEnd);
        if (!realSequence || realSequence.length === 0) {
            realSequence = "N".repeat(windowLength);
        }

        container.innerHTML = '';

        const viewerWrapper = document.createElement('div');
        viewerWrapper.className = 'franklin-viewer';

        // Biological adjacent exon resolution (strictly next in splicing sequence N+1 / N-1, NOT last exon of gene)
        const focusedExon = exons.find(ex => (targetPos >= ex.start && targetPos <= ex.end));
        const focusedIntron = introns.find(intr => (targetPos >= intr.start && targetPos <= intr.end));

        let prevExon = null;
        let nextExon = null;

        if (focusedExon) {
            prevExon = exons.find(e => e.exonNum === focusedExon.exonNum - 1) || null;
            nextExon = exons.find(e => e.exonNum === focusedExon.exonNum + 1) || null;
        } else if (focusedIntron) {
            // Between donorExon and acceptorExon
            prevExon = exons.find(e => e.exonNum === focusedIntron.donorExon) || null;
            nextExon = exons.find(e => e.exonNum === focusedIntron.acceptorExon) || null;
        } else {
            // Flanqueante: buscar exones adyacentes inmediatos en orden del transcripto
            if (!isAntisense) {
                const ups = exons.filter(e => e.end < targetPos).sort((a, b) => b.end - a.end);
                const downs = exons.filter(e => e.start > targetPos).sort((a, b) => a.start - b.start);
                prevExon = ups[0] || null;
                nextExon = downs[0] || null;
            } else {
                const ups = exons.filter(e => e.start > targetPos).sort((a, b) => a.start - b.start);
                const downs = exons.filter(e => e.end < targetPos).sort((a, b) => b.end - a.end);
                prevExon = ups[0] || null;
                nextExon = downs[0] || null;
            }
        }

        // Positions array strictly aligned in 5' -> 3' mRNA orientation
        const positions = [];
        if (!isAntisense) {
            for (let p = winStart; p <= winEnd; p++) positions.push(p);
        } else {
            // Hebra reversa: 5' está en winEnd (coordenada genómica más alta)
            for (let p = winEnd; p >= winStart; p--) positions.push(p);
        }

        // 1. Header Bar
        const headerBar = document.createElement('div');
        headerBar.className = 'franklin-header';
        headerBar.innerHTML = `
            <div class="franklin-pill">
                Chr${chromosome}:${winStart.toLocaleString()}–${winEnd.toLocaleString()} (${model.strand}) 🔍
            </div>
            <div class="franklin-nav">
                <button id="fNav-prevExon" class="btn-secondary" title="Saltar al Exón Anterior (${prevExon ? `Exón ${prevExon.exonNum}` : 'Inicio'})" ${!prevExon ? 'disabled style="opacity:0.5;"' : ''}>
                    ⏮ ${prevExon ? `Exón ${prevExon.exonNum}` : 'Inicio'}
                </button>
                <button id="fNav-prev" class="btn-secondary" title="Desplazar 10 pb hacia el 5' (río arriba)">◄ −10 pb</button>
                <button id="fNav-center" class="btn-secondary" title="Centrar en la posición de la variante">🎯 Centrar Variante</button>
                <button id="fNav-next" class="btn-secondary" title="Desplazar 10 pb hacia el 3' (río abajo)">+10 pb ►</button>
                <button id="fNav-nextExon" class="btn-secondary" title="Saltar al Próximo Exón (${nextExon ? `Exón ${nextExon.exonNum}` : 'Fin'})" ${!nextExon ? 'disabled style="opacity:0.5;"' : ''}>
                    ${nextExon ? `Exón ${nextExon.exonNum}` : 'Fin'} ⏭
                </button>
            </div>
            <div class="franklin-view-mode-toggle">
                <button id="toggleWtMutBtn" class="btn-toggle-mode ${isShowingMutantView ? 'is-mut' : 'is-wt'}">
                    ${isShowingMutantView ? '⚡ Mostrando: Variante Mutada' : '🧬 Mostrando: Secuencia Salvaje (WT)'}
                </button>
            </div>
        `;
        viewerWrapper.appendChild(headerBar);


        const viewport = document.createElement('div');
        viewport.className = 'franklin-viewport';

        // Synchronized Inner Tracks Container to guarantee 1:1 pixel alignment across all columns
        const innerTracks = document.createElement('div');
        innerTracks.className = 'f-tracks-inner';
        innerTracks.style.minWidth = `${windowLength * 30}px`;

        // 2.1 Genomic Ruler Track
        const rulerTrack = document.createElement('div');
        rulerTrack.className = 'f-ruler-track';
        rulerTrack.style.display = 'grid';
        rulerTrack.style.gridTemplateColumns = `repeat(${windowLength}, minmax(30px, 1fr))`;

        for (let i = 0; i < windowLength; i++) {
            const p = positions[i];
            const cell = document.createElement('div');
            cell.className = 'f-ruler-cell';
            if (p % 10 === 0 || i === 0 || i === windowLength - 1) {
                cell.innerHTML = `
                    <span class="f-tick">◆</span>
                    <span class="f-pos">${p.toLocaleString()}</span>
                `;
            }
            rulerTrack.appendChild(cell);
        }
        innerTracks.appendChild(rulerTrack);

        function getExonAt(pos) {
            return exons.find(e => pos >= e.start && pos <= e.end);
        }

        // 2.2 Region Header Track (Compact, non-distorting)
        const regionHeaderTrack = document.createElement('div');
        regionHeaderTrack.className = 'f-region-header-track';
        regionHeaderTrack.style.display = 'grid';
        regionHeaderTrack.style.gridTemplateColumns = `repeat(${windowLength}, minmax(30px, 1fr))`;

        for (let i = 0; i < windowLength; ) {
            const p = positions[i];
            const curEx = getExonAt(p);
            if (curEx) {
                let spanLen = 0;
                while (i + spanLen < windowLength && getExonAt(positions[i + spanLen])?.exonNum === curEx.exonNum) {
                    spanLen++;
                }
                const box = document.createElement('div');
                box.className = 'f-region-box is-exon';
                box.style.gridColumn = `span ${spanLen}`;
                box.innerHTML = `<span>[ Exón ${curEx.exonNum} ]</span>`;
                box.title = `Exón ${curEx.exonNum} (${curEx.length.toLocaleString()} pb) [Chr${chromosome}:${curEx.start.toLocaleString()}–${curEx.end.toLocaleString()}]`;
                regionHeaderTrack.appendChild(box);
                i += spanLen;
            } else {
                const curIntr = introns.find(intron => p >= intron.start && p <= intron.end);
                let spanLen = 0;
                while (i + spanLen < windowLength && !getExonAt(positions[i + spanLen])) {
                    spanLen++;
                }
                const box = document.createElement('div');
                box.className = 'f-region-box is-intron';
                box.style.gridColumn = `span ${spanLen}`;
                box.innerHTML = `<span>[ Intrón ${curIntr ? curIntr.intronNum : ''} ]</span>`;
                box.title = `Región Intrónica ${curIntr ? `[Intrón ${curIntr.intronNum}] (entre Exón ${curIntr.donorExon} y Exón ${curIntr.acceptorExon})` : ''}`;
                regionHeaderTrack.appendChild(box);
                i += spanLen;
            }
        }
        innerTracks.appendChild(regionHeaderTrack);

        // 2.3 Amino Acid & Chevron Track (Siempre en orientación biológica 5' -> 3' del ARNm)
        const aaTrack = document.createElement('div');
        aaTrack.className = 'f-aa-track';
        aaTrack.style.display = 'grid';
        aaTrack.style.gridTemplateColumns = `repeat(${windowLength}, minmax(30px, 1fr))`;

        const refLen = Math.max(1, (variant.ref || 'N').length);
        const varStart = variant.pos;
        const varEnd = variant.pos + refLen - 1;

        for (let i = 0; i < windowLength; ) {
            const p = positions[i];
            const currentExon = getExonAt(p);
            
            if (currentExon && currentExon.isCoding) {
                const offsetInExon = isAntisense ? (currentExon.end - p) : (p - currentExon.start);
                const phase = (offsetInExon + currentExon.phase) % 3;
                const remainingBasesInThisExon = isAntisense ? (p - currentExon.start + 1) : (currentExon.end - p + 1);

                if (phase === 0) {
                    if (remainingBasesInThisExon >= 3 && i + 2 < windowLength) {
                        let codonBases = "";
                        for (let c = 0; c < 3; c++) {
                            const posC = positions[i + c];
                            let baseC = (realSequence[posC - winStart] || 'N').toUpperCase();
                            if (isAntisense) baseC = reverseComplement(baseC);
                            
                            // Si la variante solapa esta base en vista mutada
                            if (isShowingMutantView && posC >= varStart && posC <= varEnd) {
                                const varOff = isAntisense ? (varEnd - posC) : (posC - varStart);
                                const senseAlt = isAntisense ? reverseComplement(variant.alt) : variant.alt;
                                baseC = senseAlt[varOff] || baseC;
                            }
                            codonBases += baseC;
                        }

                        const aa = CODON_TABLE[codonBases] || 'X';
                        const aa3 = aa1to3(aa);

                        const chevron = document.createElement('div');
                        chevron.className = `f-aa-chevron ${aa === '*' ? 'is-stop' : ''}`;
                        chevron.style.gridColumn = 'span 3';
                        chevron.innerHTML = `<span>${aa === '*' ? 'Stop' : aa3}</span>`;
                        chevron.title = `Codón ${codonBases} → ${aaFullName(aa)} (${aa3}) [5'➔3' ARNm]`;
                        aaTrack.appendChild(chevron);
                        i += 3;
                    } else if (remainingBasesInThisExon < 3) {
                        const splitCount = remainingBasesInThisExon;
                        const actualSpan = Math.min(splitCount, windowLength - i);
                        if (actualSpan <= 0) break;

                        let splitInfo = null;
                        try {
                            splitInfo = await getSplitCodonInfo(currentExon, exons, chromosome);
                        } catch (e) {
                            console.warn("Split codon lookup notice:", e);
                        }

                        const splitEl = document.createElement('div');
                        splitEl.className = 'f-aa-split';
                        splitEl.style.gridColumn = `span ${actualSpan}`;

                        if (splitInfo) {
                            splitEl.innerHTML = `<span>[${splitInfo.aa3} ${splitCount}/3 ➔]</span>`;
                            splitEl.title = `Codón dividido en límite de Exón ${currentExon.exonNum}:\n${splitInfo.tailSeq} (Exón ${splitInfo.currentExonNum}) + ${splitInfo.headSeq} (Exón ${splitInfo.nextExonNum}) = ${splitInfo.fullCodon} → ${splitInfo.name} (${splitInfo.aa3})`;
                        } else {
                            splitEl.innerHTML = `<span>[Fase ${splitCount}/3 ➔]</span>`;
                            splitEl.title = `Codón dividido en límite de Exón ${currentExon.exonNum}: ${splitCount} nt (se completa en el siguiente exón)`;
                        }

                        aaTrack.appendChild(splitEl);
                        i += actualSpan;
                    } else {
                        const partialSpan = windowLength - i;
                        const partCell = document.createElement('div');
                        partCell.className = 'f-aa-partial';
                        partCell.style.gridColumn = `span ${partialSpan}`;
                        partCell.innerHTML = `<span>...</span>`;
                        aaTrack.appendChild(partCell);
                        i += partialSpan;
                    }
                } else {
                    const midCell = document.createElement('div');
                    midCell.className = 'f-aa-mid';
                    aaTrack.appendChild(midCell);
                    i++;
                }
            } else if (!currentExon) {
                const intronCell = document.createElement('div');
                intronCell.className = 'f-intron-cell';
                intronCell.textContent = (i % 3 === 0) ? '►' : '';
                aaTrack.appendChild(intronCell);
                i++;
            } else {
                const utrCell = document.createElement('div');
                utrCell.className = 'f-utr-cell';
                aaTrack.appendChild(utrCell);
                i++;
            }
        }
        innerTracks.appendChild(aaTrack);

        // 2.4 Base Letters Track (Con soporte de hebra sense 5' -> 3')
        const baseTrack = document.createElement('div');
        baseTrack.className = 'f-base-track';
        baseTrack.style.display = 'grid';
        baseTrack.style.gridTemplateColumns = `repeat(${windowLength}, minmax(30px, 1fr))`;

        const senseRefStr = isAntisense ? reverseComplement(variant.ref) : variant.ref;
        const senseAltStr = isAntisense ? reverseComplement(variant.alt) : variant.alt;

        for (let i = 0; i < windowLength; i++) {
            const p = positions[i];
            let baseLetter = (realSequence[p - winStart] || 'N').toUpperCase();
            if (isAntisense) baseLetter = reverseComplement(baseLetter);

            const isWithinVariant = (p >= varStart && p <= varEnd);

            const baseCell = document.createElement('div');
            baseCell.className = `f-base-cell ${isWithinVariant ? 'is-variant-pos' : ''}`;

            if (isWithinVariant) {
                const offset = isAntisense ? (varEnd - p) : (p - varStart);
                if (isShowingMutantView) {
                    baseCell.classList.add('is-variant-mut-span');
                    if (senseRefStr.length > senseAltStr.length) {
                        // Deleción
                        if (offset < senseAltStr.length) {
                            const altChar = senseAltStr[offset];
                            baseCell.innerHTML = `<span class="f-variant-badge mut-badge" title="Variante mutada (${senseRefStr} > ${senseAltStr}) [Base ${offset + 1}]: ${altChar}">${altChar}</span>`;
                        } else {
                            const delChar = senseRefStr[offset] || baseLetter;
                            baseCell.innerHTML = `<span class="f-variant-badge del-badge" title="Base WT eliminada por deleción: ${delChar} en Chr${chromosome}:${p}">−</span>`;
                        }
                    } else if (senseRefStr.length < senseAltStr.length) {
                        // Inserción
                        if (offset === 0) {
                            baseCell.innerHTML = `<span class="f-variant-badge ins-badge" title="Inserción en Chr${chromosome}:${p}: ${senseRefStr} > ${senseAltStr}">+${senseAltStr.slice(senseRefStr.length)}</span>`;
                        } else {
                            baseCell.innerHTML = `<span class="f-variant-badge mut-badge" title="Inserción">+</span>`;
                        }
                    } else {
                        // Sustitución
                        const altChar = senseAltStr[offset] || senseAltStr;
                        baseCell.innerHTML = `<span class="f-variant-badge mut-badge" title="Variante Mutada: ${senseRefStr[offset] || senseRefStr} > ${altChar} en Chr${chromosome}:${p}">${altChar}</span>`;
                    }
                } else {
                    // Vista WT: Resaltar todas las bases del alelo de referencia
                    baseCell.classList.add('is-variant-wt-span');
                    const wtChar = (senseRefStr && senseRefStr !== 'N' && offset < senseRefStr.length) ? senseRefStr[offset] : baseLetter;
                    baseCell.innerHTML = `<span class="f-variant-badge wt-badge" title="Base WT implicada en la variante (${senseRefStr} > ${senseAltStr}) [Locus Chr${chromosome}:${p}, base ${offset + 1} de ${senseRefStr.length}]: ${wtChar}">${wtChar}</span>`;
                }
            } else {
                const color = BASE_COLORS[baseLetter] || BASE_COLORS['N'];
                baseCell.innerHTML = `<span class="f-base" style="color: ${color};">${baseLetter}</span>`;
                baseCell.title = `Chr${chromosome}:${p} = ${baseLetter}${isAntisense ? ' (Sense 5\'➔3\')' : ''}`;
            }

            baseTrack.appendChild(baseCell);
        }
        innerTracks.appendChild(baseTrack);

        // 2.5 Splice Junction Demarcation Lines
        const junctionsInWindow = [];
        exons.forEach(ex => {
            if (ex.end >= winStart && ex.end <= winEnd) {
                junctionsInWindow.push({ pos: ex.end, type: isAntisense ? 'acceptor' : 'donor', exonNum: ex.exonNum });
            }
            if (ex.start >= winStart && ex.start <= winEnd) {
                junctionsInWindow.push({ pos: ex.start, type: isAntisense ? 'donor' : 'acceptor', exonNum: ex.exonNum });
            }
        });

        junctionsInWindow.forEach(j => {
            let jCol = -1;
            if (!isAntisense) {
                jCol = j.pos - winStart;
            } else {
                jCol = winEnd - j.pos;
            }
            if (jCol >= 0 && jCol < windowLength) {
                const line = document.createElement('div');
                line.className = 'f-junction-line';
                line.style.left = `calc(${((jCol + 1) / windowLength) * 100}% - 1px)`;
                line.title = j.type === 'donor' 
                    ? `Límite Exón ${j.exonNum} / Intrón (Dador)`
                    : `Límite Intrón / Exón ${j.exonNum} (Aceptor)`;
                innerTracks.appendChild(line);
            }
        });

        // 2.6 Variant Marker Line / Area (spans all affected reference positions)
        if (varEnd >= winStart && varStart <= winEnd) {
            let startCol, endCol;
            if (!isAntisense) {
                startCol = Math.max(0, varStart - winStart);
                endCol = Math.min(windowLength - 1, varEnd - winStart);
            } else {
                startCol = Math.max(0, winEnd - varEnd);
                endCol = Math.min(windowLength - 1, winEnd - varStart);
            }
            const spanBases = endCol - startCol + 1;
            const vLine = document.createElement('div');
            vLine.className = 'f-variant-marker-line';
            vLine.style.left = `calc(${((startCol) / windowLength) * 100}%)`;
            vLine.style.width = `calc(${((spanBases) / windowLength) * 100}%)`;
            vLine.title = `Locus de la variante: Chr${chromosome}:${varStart}${refLen > 1 ? `–${varEnd}` : ''} (${variant.ref} > ${variant.alt})`;
            innerTracks.appendChild(vLine);
        }

        viewport.appendChild(innerTracks);
        viewerWrapper.appendChild(viewport);
        container.appendChild(viewerWrapper);

        // Event Listeners for Nav buttons
        const navPrevExon = viewerWrapper.querySelector('#fNav-prevExon');
        const navPrev = viewerWrapper.querySelector('#fNav-prev');
        const navCenter = viewerWrapper.querySelector('#fNav-center');
        const navNext = viewerWrapper.querySelector('#fNav-next');
        const navNextExon = viewerWrapper.querySelector('#fNav-nextExon');
        const toggleWtMutBtn = viewerWrapper.querySelector('#toggleWtMutBtn');

        if (navPrevExon && prevExon) {
            navPrevExon.addEventListener('click', (e) => {
                e.stopPropagation();
                const jumpPos = isAntisense ? prevExon.end : prevExon.start;
                renderFranklinBaseViewer(container, model, jumpPos);
            });
        }

        if (navPrev) {
            navPrev.addEventListener('click', (e) => {
                e.stopPropagation();
                // Desplazar 10 pb hacia el 5' (río arriba)
                const delta = isAntisense ? 10 : -10;
                renderFranklinBaseViewer(container, model, Math.max(1, targetPos + delta));
            });
        }

        if (navCenter) {
            navCenter.addEventListener('click', (e) => {
                e.stopPropagation();
                renderFranklinBaseViewer(container, model, variant.pos);
            });
        }

        if (navNext) {
            navNext.addEventListener('click', (e) => {
                e.stopPropagation();
                // Desplazar 10 pb hacia el 3' (río abajo)
                const delta = isAntisense ? -10 : 10;
                renderFranklinBaseViewer(container, model, targetPos + delta);
            });
        }

        if (navNextExon && nextExon) {
            navNextExon.addEventListener('click', (e) => {
                e.stopPropagation();
                const jumpPos = isAntisense ? nextExon.end : nextExon.start;
                renderFranklinBaseViewer(container, model, jumpPos);
            });
        }

        if (toggleWtMutBtn) {
            toggleWtMutBtn.addEventListener('click', (e) => {
                e.stopPropagation();
                renderFranklinBaseViewer(container, model, targetPos, !isShowingMutantView);
            });
        }
    } catch (err) {
        console.error("Error en renderFranklinBaseViewer:", err);
        container.innerHTML = `
            <div class="nmd-decision-box trigger-nmd" style="margin: 12px 0;">
                <div class="nmd-title"><span>⚠️ No se pudo renderizar el visor micro de bases:</span></div>
                <p class="nmd-description">${err.message || 'Error inesperado al generar las pistas de secuencia.'}</p>
            </div>
        `;
    }
}

/**
 * ═══════════════════════════════════════════════════════════════
 * SYNCHRONIZED GRID TRACK RENDERER — SEGMENT-AWARE
 *
 * Accepts an array of segments. Each segment:
 *   { seq, label, type: 'exon'|'intron'|'deleted', startPhase, exonNum }
 * Options:
 *   - showAminoAcids: boolean (false for Pre-mRNA nucleotide-only stages)
 *   - stopAtCodonStop: boolean (true to truncate strictly at the first STOP codon)
 * ═══════════════════════════════════════════════════════════════
 */
export function renderSynchronizedGridTrack(segments = [], {
    showAminoAcids = true,
    stopAtCodonStop = false,
    isFrameshift = false,
    frameshiftStartIndex = null,
    variantPosIndex = null,
    variantRef = 'N',
    variantAlt = 'N',
    junctionIndexes = []
} = {}) {

    // Flatten segments into a per-character metadata array
    let chars = [];
    for (let sIdx = 0; sIdx < segments.length; sIdx++) {
        const seg = segments[sIdx];
        const seq = (seg.seq || '').toUpperCase();
        for (let i = 0; i < seq.length; i++) {
            chars.push({
                base: seq[i] || 'N',
                type: seg.type || 'exon',
                segLabel: seg.label || '',
                sIdx,
                localIdx: i
            });
        }
    }

    if (chars.length === 0) return document.createElement('div');

    const viewport = document.createElement('div');
    viewport.className = 'franklin-viewport sync-grid-viewport';
    viewport.style.position = 'relative';

    let gIdx = 0; // global index into chars[]
    let rollingPhase = 0;
    let stopCutoffIndex = null;

    // === 1. AMINO ACID TRACK (Only if showAminoAcids is true) ===
    let aaTrack = null;
    if (showAminoAcids) {
        aaTrack = document.createElement('div');
        aaTrack.className = 'f-aa-track sync-grid-row';
        aaTrack.style.display = 'grid';

        for (let sIdx = 0; sIdx < segments.length; sIdx++) {
            if (stopCutoffIndex !== null) break;

            const seg = segments[sIdx];
            const seq = (seg.seq || '').toUpperCase();
            const segLen = seq.length;
            if (segLen === 0) { continue; }

            const phase = (seg.startPhase !== undefined && seg.startPhase !== null)
                ? seg.startPhase
                : rollingPhase;

            const isIntronSeg = (seg.type === 'intron');
            const isDeletedSeg = (seg.type === 'deleted');

            // --- DELETED SEGMENT (e.g. omitted exon in skipping) ---
            if (isDeletedSeg) {
                const delCell = document.createElement('div');
                delCell.className = 'f-aa-deleted';
                delCell.style.gridColumn = `span ${segLen}`;
                delCell.style.display = 'flex';
                delCell.style.alignItems = 'center';
                delCell.style.justifyContent = 'center';
                delCell.style.fontSize = '0.7rem';
                delCell.style.color = 'var(--accent-rose)';
                delCell.style.textDecoration = 'line-through';
                delCell.style.opacity = '0.65';
                delCell.innerHTML = `<span>${seg.label || 'Omitido'}</span>`;
                aaTrack.appendChild(delCell);
                gIdx += segLen;
                rollingPhase = (phase + segLen) % 3;
                continue;
            }

            // --- INTRON SEGMENT ---
            if (isIntronSeg) {
                let ip = 0;
                const leadBasesCount = phase > 0 ? ((3 - phase) % 3) : 0;
                if (leadBasesCount > 0 && ip < segLen) {
                    const leadCell = document.createElement('div');
                    leadCell.className = 'f-aa-mid f-intron-aa-mid';
                    leadCell.style.gridColumn = `span ${leadBasesCount}`;
                    leadCell.style.background = 'rgba(168,85,247,0.1)';
                    aaTrack.appendChild(leadCell);
                    ip += leadBasesCount;
                }
                while (ip + 2 < segLen) {
                    const codon = seq.substring(ip, ip + 3);
                    const aa = CODON_TABLE[codon] || 'X';
                    const aa3 = aa1to3(aa);
                    const isStop = (aa === '*');
                    const chevron = document.createElement('div');
                    chevron.style.gridColumn = 'span 3';
                    if (isStop) {
                        chevron.className = 'f-aa-chevron is-stop';
                        chevron.style.cssText += ';background:rgba(239,68,68,0.45);border-color:var(--accent-rose);color:#fff;box-shadow:0 0 12px rgba(239,68,68,0.9);';
                        chevron.innerHTML = `<span style="font-weight:900;">🛑 STOP *</span>`;
                        chevron.title = `Codón Stop Prematuro en intrón retenido: ${codon}`;
                        aaTrack.appendChild(chevron);
                        ip += 3;

                        // Strict termination at STOP
                        if (stopAtCodonStop) {
                            stopCutoffIndex = gIdx + ip;
                            break;
                        }
                    } else {
                        chevron.className = 'f-aa-chevron f-intron-aa-chevron';
                        chevron.style.cssText += ';background:rgba(168,85,247,0.2);border-color:#a855f7;color:#c084fc;';
                        chevron.innerHTML = `<span>${aa3}</span>`;
                        chevron.title = `Codón intrónico: ${codon} → ${aaFullName(aa)} (${aa3})`;
                        aaTrack.appendChild(chevron);
                        ip += 3;
                    }
                }
                if (stopCutoffIndex === null && ip < segLen) {
                    const rem = document.createElement('div');
                    rem.className = 'f-aa-partial f-intron-aa-mid';
                    rem.style.gridColumn = `span ${segLen - ip}`;
                    rem.style.background = 'rgba(168,85,247,0.08)';
                    rem.innerHTML = `<span style="opacity:0.4">···</span>`;
                    aaTrack.appendChild(rem);
                }
                gIdx += (stopCutoffIndex !== null ? ip : segLen);
                rollingPhase = (phase + segLen) % 3;
                continue;
            }

            // --- EXON SEGMENT ---
            let ep = 0;

            // 1. If segment begins with a split-codon completion from previous exon (phase = 1 or 2)
            const leadBasesCount = phase > 0 ? ((3 - phase) % 3) : 0;
            if (leadBasesCount > 0 && ep < segLen) {
                const completingBases = seq.substring(0, leadBasesCount);
                const prevSeg = sIdx > 0 ? segments[sIdx - 1] : null;
                const prevSeq = prevSeg ? (prevSeg.seq || '').toUpperCase() : '';
                const prevTail = prevSeq.substring(Math.max(0, prevSeq.length - phase));
                const fullCodon = (prevTail + completingBases).toUpperCase();
                const aa = CODON_TABLE[fullCodon] || 'X';
                const aa3 = aa1to3(aa);

                const splitEnd = document.createElement('div');
                splitEnd.className = 'f-aa-split';
                splitEnd.style.gridColumn = `span ${leadBasesCount}`;
                splitEnd.innerHTML = `<span>[➔ ${leadBasesCount}/3 ${aa3}]</span>`;
                splitEnd.title = `Codón dividido completado: ${prevTail || '...'}+${completingBases} = ${fullCodon} → ${aaFullName(aa)} (${aa3})`;
                aaTrack.appendChild(splitEnd);
                ep += leadBasesCount;
            }

            // 2. Translate full codons
            while (ep + 2 < segLen) {
                const remaining = segLen - ep;
                if (remaining < 3) break; // handle below as trailing split

                const codon = seq.substring(ep, ep + 3);
                const aa = CODON_TABLE[codon] || 'X';
                const aa3 = aa1to3(aa);
                const isStop = (aa === '*');
                const globalCodonStart = gIdx + ep;
                const isAberrant = isFrameshift && frameshiftStartIndex !== null && globalCodonStart >= frameshiftStartIndex;

                const chevron = document.createElement('div');
                chevron.style.gridColumn = 'span 3';

                if (isStop) {
                    chevron.className = 'f-aa-chevron is-stop';
                    chevron.style.cssText += ';background:rgba(239,68,68,0.45);border-color:var(--accent-rose);color:#fff;box-shadow:0 0 12px rgba(239,68,68,0.9);';
                    chevron.innerHTML = `<span style="font-weight:900;">🛑 STOP *</span>`;
                    chevron.title = `Codón Stop Prematuro: ${codon}`;
                    aaTrack.appendChild(chevron);
                    ep += 3;

                    // Strict termination at STOP
                    if (stopAtCodonStop) {
                        stopCutoffIndex = gIdx + ep;
                        break;
                    }
                } else if (isAberrant) {
                    chevron.className = 'f-aa-chevron';
                    chevron.style.cssText += ';background:rgba(245,158,11,0.25);border-color:var(--accent-amber);color:#fbbf24;';
                    chevron.innerHTML = `<span>${aa3}</span>`;
                    chevron.title = `Marco aberrante: ${codon} → ${aaFullName(aa)} (${aa3})`;
                    aaTrack.appendChild(chevron);
                    ep += 3;
                } else {
                    chevron.className = 'f-aa-chevron';
                    chevron.innerHTML = `<span>${aa3}</span>`;
                    chevron.title = `${seg.label ? seg.label + ': ' : ''}${codon} → ${aaFullName(aa)} (${aa3})`;
                    aaTrack.appendChild(chevron);
                    ep += 3;
                }
            }

            if (stopCutoffIndex !== null) {
                gIdx += ep;
                break;
            }

            // 3. Trailing split codon at end of this exon segment
            const endSplitCount = (phase + segLen) % 3;
            if (endSplitCount > 0 && ep < segLen) {
                const tailBases = seq.substring(ep);
                const tailLen = tailBases.length;
                const nextSeg = sIdx < segments.length - 1 ? segments[sIdx + 1] : null;
                const nextSeq = nextSeg ? (nextSeg.seq || '').toUpperCase() : '';
                const headBases = nextSeq.substring(0, 3 - tailLen);
                const fullCodon = (tailBases + headBases).toUpperCase();
                const aa = CODON_TABLE[fullCodon] || 'X';
                const aa3 = aa1to3(aa);

                const splitStart = document.createElement('div');
                splitStart.className = 'f-aa-split';
                splitStart.style.gridColumn = `span ${tailLen}`;
                splitStart.innerHTML = `<span>[${aa3} ${tailLen}/3 ➔]</span>`;
                splitStart.title = `Codón dividido: ${tailBases} (${tailLen}nt aquí) + ${headBases || '...'} = ${fullCodon} → ${aaFullName(aa)} (${aa3})`;
                aaTrack.appendChild(splitStart);
                ep += tailLen;
            } else if (ep < segLen) {
                const rem = document.createElement('div');
                rem.className = 'f-aa-partial';
                rem.style.gridColumn = `span ${segLen - ep}`;
                rem.innerHTML = `<span>···</span>`;
                aaTrack.appendChild(rem);
            }

            gIdx += segLen;
            rollingPhase = (phase + segLen) % 3;
        }
    }

    // Apply strict STOP cutoff to character array if triggered
    if (stopCutoffIndex !== null && stopCutoffIndex < chars.length) {
        chars = chars.slice(0, stopCutoffIndex);
    }

    const finalLen = chars.length;
    if (aaTrack) {
        aaTrack.style.gridTemplateColumns = `repeat(${finalLen}, 32px)`;
        viewport.appendChild(aaTrack);
    }

    // === 2. BASE LETTERS TRACK ===
    const baseTrack = document.createElement('div');
    baseTrack.className = 'f-base-track sync-grid-row';
    baseTrack.style.display = 'grid';
    baseTrack.style.gridTemplateColumns = `repeat(${finalLen}, 32px)`;

    for (let i = 0; i < finalLen; i++) {
        const ch = chars[i];
        const b = ch.base;
        const color = BASE_COLORS[b] || '#fff';
        const isVar = (i === variantPosIndex);
        const isIntron = (ch.type === 'intron');
        const isDeleted = (ch.type === 'deleted');

        const baseCell = document.createElement('div');
        baseCell.className = [
            'f-base-cell',
            isVar ? 'is-variant-pos' : '',
            isIntron ? 'is-intron-base' : '',
            isDeleted ? 'is-deleted-base' : ''
        ].join(' ').trim();

        if (isVar) {
            baseCell.innerHTML = `<span class="f-variant-badge mut-badge" title="Variante Mutada: ${variantRef}>${variantAlt}">${variantAlt || b}</span>`;
        } else {
            const regionLabel = isIntron ? ' (Intrón)' : isDeleted ? ' (Omitido)' : ` (${ch.segLabel || 'Exón'})`;
            baseCell.innerHTML = `<span class="f-base" style="color: ${color};">${b}</span>`;
            baseCell.title = `Pos +${i + 1}: ${b}${regionLabel}`;
        }
        baseTrack.appendChild(baseCell);
    }
    viewport.appendChild(baseTrack);

    // === 3. VERTICAL DEMARCATION LINES ===
    for (const jIdx of junctionIndexes) {
        if (jIdx > 0 && jIdx <= finalLen) {
            const linePct = (jIdx / finalLen) * 100;
            const jLine = document.createElement('div');
            jLine.className = 'f-junction-line';
            jLine.style.left = `calc(${linePct}% - 1px)`;
            jLine.style.background = '#38bdf8';
            jLine.style.boxShadow = '0 0 8px rgba(56,189,248,0.8)';
            viewport.appendChild(jLine);
        }
    }

    return viewport;
}

/**
 * ═══════════════════════════════════════════════════════════════
 * VISOR DUAL DIDÁCTICO PASO A PASO (CANÓNICO WT vs MODIFICADO)
 * ═══════════════════════════════════════════════════════════════
 */
export function renderComparisonSplicingViewer(containerEl, {
    mode = 'skipping',
    prevExon = null,
    targetExon = null,
    nextExon = null,
    skippedLen = 0,
    isFrameshift = false,
    frameshiftShift = 0,
    seqPrev = "CAGTACCAGTTGAC",
    seqPrevPhase = 0,
    seqTarget = "TTTGAAAGTGATGAA",
    seqTargetTail = null,
    seqTargetPhase = 0,
    seqNext = "TTAGCTGAATTGGAC",
    seqNextPhase = 0,
    crypticDelta = 4,
    isIntronCryptic = true,
    crypticIntronSeq = null,
    variantPosIndex = null,
    retainedIntronSeq = "GTAAGTTAGCTAATGACTTGACCA",
    variant = { ref: 'G', alt: 'A', pos: 0 }
}) {
    if (!containerEl) return;
    containerEl.style.display = 'block';
    containerEl.innerHTML = '';

    const wrap = document.createElement('div');
    wrap.className = 'splicing-comparison-container';

    // Helper: build a named box with header + track
    function makeNamedBox(cls, headerHtml, trackEl) {
        const box = document.createElement('div');
        box.className = `track-named-box ${cls}`;
        const hdr = document.createElement('div');
        hdr.className = 'track-box-header';
        hdr.innerHTML = headerHtml;
        box.appendChild(hdr);
        box.appendChild(trackEl);
        return box;
    }

    // Helper: make a stage section card
    function makeStageCard(cls, badgeHtml, subtextHtml) {
        const sec = document.createElement('div');
        sec.className = `comparison-stage-card ${cls}`;
        sec.innerHTML = `
            <div class="stage-title">
                ${badgeHtml}
                <span class="stage-subtext">${subtextHtml}</span>
            </div>
        `;
        return sec;
    }

    // ─────────────────────────────────────────────────────────────
    // CASO 1: SALTO DE EXÓN COMPLETO (TAREA 2)
    // Estructura directa:
    // [Final Exón anterior] (secuencia)
    // "Se perdería:"
    // [Inicio Exón saltado] --//-- [Final Exón saltado]
    // [Inicio Exón siguiente]
    // Resultado: [Empalme final Exón anterior unido a Exón siguiente]
    // ─────────────────────────────────────────────────────────────
    if (mode === 'skipping') {
        const pNum = prevExon ? prevExon.exonNum : 1;
        const tNum = targetExon ? targetExon.exonNum : 2;
        const nNum = nextExon ? nextExon.exonNum : 3;

        const prevPhase   = seqPrevPhase   !== null ? seqPrevPhase   : 0;
        const targetPhase = seqTargetPhase !== null ? seqTargetPhase : (targetExon?.phase ?? 0);
        const nextPhase   = seqNextPhase   !== null ? seqNextPhase   : (nextExon?.phase ?? 0);

        const wtSection = makeStageCard('wt-stage',
            `<span class="badge" style="background:rgba(16,185,129,0.18);border:1px solid #10b981;color:#34d399;">🧬 1. Estructura Canónica Salvaje (WT)</span>`,
            `Límites de los exones originales y región que se omite`
        );
        const wtFlex = document.createElement('div');
        wtFlex.className = 'sync-tracks-column';

        // 1. [Final Exón anterior]
        const ex1Track = renderSynchronizedGridTrack(
            [{ seq: seqPrev, label: `Exón ${pNum}`, type: 'exon', startPhase: prevPhase }],
            { junctionIndexes: [] }
        );
        wtFlex.appendChild(makeNamedBox('',
            `<span>[ Final Exón <strong>${pNum}</strong> ]</span><span class="phase-tag">Salida: Fase ${prevExon?.endPhase ?? 0}</span>`,
            ex1Track
        ));

        // 2. Texto: "Se perdería:"
        const lossDivider = document.createElement('div');
        lossDivider.className = 'skipping-loss-divider';
        lossDivider.innerHTML = `<span>⚠️ Se perdería en el transcripto maduro:</span>`;
        wtFlex.appendChild(lossDivider);

        // 3. [Inicio Exón saltado] --//-- [Final Exón saltado] (con elisor si > 30 pb)
        if (skippedLen > 30 && seqTargetTail) {
            const exTHeadTrack = renderSynchronizedGridTrack(
                [{ seq: seqTarget, label: `Exón ${tNum} (Inicio)`, type: 'exon', startPhase: targetPhase }],
                { junctionIndexes: [] }
            );
            wtFlex.appendChild(makeNamedBox('is-skipped-box',
                `<span style="color:var(--accent-rose);">[ Inicio Exón <strong>${tNum}</strong> (15 pb) ] — Omitido</span>
                 <span class="phase-tag">Entrada: Fase ${targetPhase}</span>`,
                exTHeadTrack
            ));

            const elisorEl = document.createElement('div');
            elisorEl.className = 'exon-elisor-separator';
            elisorEl.innerHTML = `<span class="elisor-line">──//──</span> <span class="elisor-badge">... // Exón ${tNum} Omitido (${skippedLen} pb en total) // ...</span> <span class="elisor-line">──//──</span>`;
            wtFlex.appendChild(elisorEl);

            const exTTailTrack = renderSynchronizedGridTrack(
                [{ seq: seqTargetTail, label: `Exón ${tNum} (Final)`, type: 'exon', startPhase: 0 }],
                { junctionIndexes: [] }
            );
            wtFlex.appendChild(makeNamedBox('is-skipped-box',
                `<span style="color:var(--accent-rose);">[ Final Exón <strong>${tNum}</strong> (15 pb) ] — Omitido</span>
                 <span class="phase-tag">Salida: Fase ${targetExon?.endPhase ?? 0}</span>`,
                exTTailTrack
            ));
        } else {
            const exTTrack = renderSynchronizedGridTrack(
                [{ seq: seqTarget, label: `Exón ${tNum}`, type: 'exon', startPhase: targetPhase }],
                { junctionIndexes: [] }
            );
            wtFlex.appendChild(makeNamedBox('is-skipped-box',
                `<span style="color:var(--accent-rose);">[ Exón <strong>${tNum}</strong> (${skippedLen} pb) ] — Omitido</span>
                 <span class="badge badge-danger" style="font-size:0.65rem;">Se omite</span>`,
                exTTrack
            ));
        }

        const nextDivider = document.createElement('div');
        nextDivider.className = 'intron-flow-separator';
        nextDivider.innerHTML = `<span>──►</span><span class="intron-split-badge">Intrón (-//-)</span><span>──►</span>`;
        wtFlex.appendChild(nextDivider);

        // 4. [Inicio Exón siguiente]
        const ex2Track = renderSynchronizedGridTrack(
            [{ seq: seqNext, label: `Exón ${nNum}`, type: 'exon', startPhase: nextPhase }],
            { junctionIndexes: [] }
        );
        wtFlex.appendChild(makeNamedBox('',
            `<span>[ Inicio Exón <strong>${nNum}</strong> ]</span><span class="phase-tag">Entrada: Fase ${nextPhase}</span>`,
            ex2Track
        ));

        wtSection.appendChild(wtFlex);
        wrap.appendChild(wtSection);

        // 5. [Resultado: Empalme final Exón anterior unido al Exón siguiente]
        const mutSection = makeStageCard('mut-stage',
            `<span class="badge ${isFrameshift ? 'badge-danger' : 'badge-warning'}">✂ 2. Resultado: ARNm Empalmado Final</span>`,
            `Empalme directo: <strong>[ Exón ${pNum} ] ✂ [ Exón ${nNum} ]</strong> (${isFrameshift ? `Frameshift +${frameshiftShift} nt → Stop Prematuro 🛑` : '✅ In-Frame: 3n conservado'})`
        );
        const mutFlex = document.createElement('div');
        mutFlex.className = 'sync-tracks-column';

        const joinedTrack = renderSynchronizedGridTrack([
            { seq: seqPrev, label: `Exón ${pNum}`, type: 'exon', startPhase: prevPhase },
            { seq: seqNext, label: `Exón ${nNum}`, type: 'exon', startPhase: null }
        ], {
            junctionIndexes: [seqPrev.length],
            isFrameshift: isFrameshift,
            frameshiftStartIndex: isFrameshift ? seqPrev.length : null
        });
        mutFlex.appendChild(makeNamedBox('full-joined-box',
            `<span>[ Resultado: Exón <strong>${pNum}</strong> ✂ Exón <strong>${nNum}</strong> ] — Unión final sin Exón ${tNum}</span>
             <span class="badge ${isFrameshift ? 'badge-danger' : 'badge-neutral'}">${isFrameshift ? '⚠️ Frameshift' : '✅ In-Frame'}</span>`,
            joinedTrack
        ));

        mutSection.appendChild(mutFlex);
        wrap.appendChild(mutSection);
    }

    // ─────────────────────────────────────────────────────────────
    // CASO 2: ACTIVACIÓN DE SITIOS CRÍPTICOS (TAREA 3)
    // Paso 1: Estructura Primaria (Pre-ARNm con región intrónica sombreada)
    // Paso 2: Empalme Mutado (ARNm unido con bloques diferenciados y traducción)
    // ─────────────────────────────────────────────────────────────
    else if (mode === 'cryptic') {
        const dNum = prevExon ? prevExon.exonNum : 1;
        const aNum = nextExon ? nextExon.exonNum : 2;
        const donorPhase = seqPrevPhase !== null ? seqPrevPhase : 0;

        if (isIntronCryptic) {
            const crypticIntronChunk = crypticIntronSeq || ("GTAAGTTCCAGTGAC").substring(0, crypticDelta);
            const newSiteIdx = seqPrev.length + crypticDelta;

            // ── PASO 1: ESTRUCTURA PRIMARIA (PRE-ARNM / NUCLEÓTIDOS CON INTRÓN SOMBREADO) ──
            const step1Card = makeStageCard('wt-stage',
                `<span class="badge" style="background:rgba(148,163,184,0.18);border:1px solid #64748b;color:#cbd5e1;">🧬 1. Estructura Primaria (Genómica / Pre-ARNm)</span>`,
                `Límite canónico original inactivado y fragmento intrónico que se incorporará (+${crypticDelta} pb)`
            );
            const rawTrack = renderSynchronizedGridTrack([
                { seq: seqPrev,            label: `Exón ${dNum} Dador`,              type: 'exon',   startPhase: donorPhase },
                { seq: crypticIntronChunk, label: `Intrón (+${crypticDelta} pb)`,   type: 'intron', startPhase: null }
            ], {
                showAminoAcids: false,
                junctionIndexes: [seqPrev.length, newSiteIdx],
                variantPosIndex: variantPosIndex,
                variantRef: variant.ref,
                variantAlt: variant.alt
            });
            step1Card.appendChild(makeNamedBox('',
                `<span>[ Exón <strong>${dNum}</strong> Dador ] ──► <span style="color:#c084fc;">[ +${crypticDelta} pb Intrónicas retenidas ]</span> ──► ⚡ Sitio Críptico</span>
                 <div style="display:flex;gap:6px;">
                   <span class="site-tag original-site-tag">📍 Límite Original</span>
                   <span class="site-tag new-site-tag">⚡ Nuevo Sitio Críptico</span>
                   <span class="pre-mrna-badge">Pre-Traducción</span>
                 </div>`,
                rawTrack
            ));
            wrap.appendChild(step1Card);

            // ── PASO 2: EMPALME MUTADO (ARNM ENSAMBLADO Y TRADUCIDO) ──
            const step2Card = makeStageCard('mut-stage',
                `<span class="badge ${isFrameshift ? 'badge-danger' : 'badge-warning'}">✂ 2. Empalme Mutado y Traducción del Marco Reconstruido</span>`,
                `Secuencia ya unida: [Exón ${dNum}] + [+${crypticDelta} pb Intrón] + [Exón ${aNum}]. Traducción continua y nuevos aminoácidos (${isFrameshift ? `Frameshift +${frameshiftShift} nt` : '✅ In-Frame'})`
            );
            const splicedTrack = renderSynchronizedGridTrack([
                { seq: seqPrev,            label: `Exón ${dNum} Dador`,              type: 'exon',   startPhase: donorPhase },
                { seq: crypticIntronChunk, label: `Intrón +${crypticDelta}pb`,       type: 'intron', startPhase: null },
                { seq: seqNext,            label: `Exón ${aNum} Aceptor`,            type: 'exon',   startPhase: null }
            ], {
                showAminoAcids: true,
                junctionIndexes: [seqPrev.length, newSiteIdx],
                variantPosIndex: variantPosIndex,
                variantRef: variant.ref,
                variantAlt: variant.alt,
                isFrameshift: isFrameshift,
                frameshiftStartIndex: isFrameshift ? seqPrev.length : null
            });
            step2Card.appendChild(makeNamedBox('is-cryptic-box',
                `<span>[ Exón <strong>${dNum}</strong> ] + <span style="color:var(--accent-amber);">[ +${crypticDelta} pb Intrón ]</span> ✂ [ Exón <strong>${aNum}</strong> ]</span>
                 <span class="badge ${isFrameshift ? 'badge-danger' : 'badge-neutral'}">${isFrameshift ? '⚠️ Frameshift Activo' : '✅ In-Frame'}</span>`,
                splicedTrack
            ));
            wrap.appendChild(step2Card);

        } else {
            // Exon Cryptic Deletion (-Delta pb)
            const keptLen = Math.max(3, seqPrev.length - crypticDelta);
            const keptSeq = seqPrev.substring(0, keptLen);

            // Paso 1
            const step1Card = makeStageCard('wt-stage',
                `<span class="badge" style="background:rgba(148,163,184,0.18);border:1px solid #64748b;color:#cbd5e1;">🧬 1. Estructura Primaria del Exón Dador</span>`,
                `Sitio críptico interno: ${crypticDelta} pb exónicas serán excluidas del ARNm maduro`
            );
            const rawTrack = renderSynchronizedGridTrack([
                { seq: keptSeq, label: `Exón ${dNum} Conservado`, type: 'exon', startPhase: donorPhase },
                { seq: seqPrev.substring(keptLen), label: `Eliminado (-${crypticDelta} pb)`, type: 'deleted', startPhase: null }
            ], {
                showAminoAcids: false,
                junctionIndexes: [keptLen],
                variantPosIndex: variantPosIndex,
                variantRef: variant.ref,
                variantAlt: variant.alt
            });
            step1Card.appendChild(makeNamedBox('',
                `<span>[ Exón <strong>${dNum}</strong> ] ──► ✂ Sitio Críptico Interno (−${crypticDelta} pb)</span>
                 <span class="site-tag new-site-tag">✂ Corte Interno</span>`,
                rawTrack
            ));
            wrap.appendChild(step1Card);

            // Paso 2
            const step2Card = makeStageCard('mut-stage',
                `<span class="badge ${isFrameshift ? 'badge-danger' : 'badge-warning'}">✂ 2. Empalme Mutado Acortado y Traducción</span>`,
                `[Exón ${dNum} acortado] unido a [Exón ${aNum}]. Traducción tras la deleción (${isFrameshift ? `Frameshift +${frameshiftShift} nt` : '✅ In-Frame'})`
            );
            const splicedTrack = renderSynchronizedGridTrack([
                { seq: keptSeq, label: `Exón ${dNum} Acortado`, type: 'exon', startPhase: donorPhase },
                { seq: seqNext, label: `Exón ${aNum} Aceptor`,   type: 'exon', startPhase: null }
            ], {
                showAminoAcids: true,
                junctionIndexes: [keptLen],
                variantPosIndex: variantPosIndex,
                variantRef: variant.ref,
                variantAlt: variant.alt,
                isFrameshift: isFrameshift,
                frameshiftStartIndex: isFrameshift ? keptLen : null
            });
            step2Card.appendChild(makeNamedBox('is-cryptic-box',
                `<span>[ Exón <strong>${dNum}</strong> (−${crypticDelta} pb) ] ✂ [ Exón <strong>${aNum}</strong> ]</span>
                 <span class="badge ${isFrameshift ? 'badge-danger' : 'badge-neutral'}">${isFrameshift ? '⚠️ Frameshift' : '✅ In-Frame'}</span>`,
                splicedTrack
            ));
            wrap.appendChild(step2Card);
        }
    }

    // ─────────────────────────────────────────────────────────────
    // CASO 3: RETENCIÓN DE INTRÓN COMPLETO (TAREA 4)
    // Marco de lectura continuo heredado del exón dador, corte estricto en STOP
    // ─────────────────────────────────────────────────────────────
    else if (mode === 'retention') {
        const dNum = prevExon ? prevExon.exonNum : 1;
        const aNum = nextExon ? nextExon.exonNum : 2;

        const donorPhase = seqPrevPhase !== null ? seqPrevPhase : 0;
        const acceptorPhase = seqNextPhase !== null ? seqNextPhase : (nextExon?.phase ?? 0);

        // === WT Stage: normal canonical splicing ===
        const wtSection = makeStageCard('wt-stage',
            `<span class="badge" style="background:rgba(16,185,129,0.18);border:1px solid #10b981;color:#34d399;">🧬 1. Splicing Canónico Salvaje (WT)</span>`,
            `El intrón es removido normalmente y los exones se unen directamente conservando el marco`
        );

        const wtTrack = renderSynchronizedGridTrack([
            { seq: seqPrev, label: `Exón ${dNum} Dador`,   type: 'exon', startPhase: donorPhase },
            { seq: seqNext, label: `Exón ${aNum} Aceptor`, type: 'exon', startPhase: acceptorPhase }
        ], {
            showAminoAcids: true,
            junctionIndexes: [seqPrev.length]
        });
        wtSection.appendChild(makeNamedBox('',
            `<span>[ Exón <strong>${dNum}</strong> ] ✂ [ Exón <strong>${aNum}</strong> ] — Empalme canónico normal</span>
             <span class="phase-tag">Proteína Canónica Intacta</span>`,
            wtTrack
        ));
        wrap.appendChild(wtSection);

        // === Mutant Stage: exon + retained intron with inherited continuous frame and strict STOP cutoff ===
        let retSeq = retainedIntronSeq || "GTAAGTTAGCTAATGACTTGACCA";
        const hasStop = /TGA|TAA|TAG/.test(retSeq.toUpperCase());
        if (!hasStop) {
            retSeq = retSeq.substring(0, 9) + "TGAACT" + retSeq.substring(9);
        }

        const mutSection = makeStageCard('mut-stage',
            `<span class="badge badge-danger">🛑 2. Retención de Intrón — Traducción Continua hasta Codón STOP (PTC)</span>`,
            `El marco de lectura del Exón ${dNum} continúa hacia el intrón retenido y se detiene estrictamente en el primer <strong>🛑 STOP *</strong>`
        );

        const mutTrack = renderSynchronizedGridTrack([
            { seq: seqPrev, label: `Exón ${dNum} Dador`,  type: 'exon',   startPhase: donorPhase },
            { seq: retSeq,  label: `Intrón Retenido`,      type: 'intron', startPhase: null }
        ], {
            showAminoAcids: true,
            stopAtCodonStop: true, // Truncates strictly at STOP
            junctionIndexes: [seqPrev.length],
            variantPosIndex: variantPosIndex !== null ? variantPosIndex : seqPrev.length,
            variantRef: variant.ref || 'G',
            variantAlt: variant.alt || 'A'
        });
        mutSection.appendChild(makeNamedBox('is-retained-box',
            `<span>[ Exón <strong>${dNum}</strong> ] ──► <span style="color:#c084fc;">[ Intrón Retenido ]</span> ──► <span style="color:var(--accent-rose);font-weight:700;">🛑 STOP * (Fin de Traducción)</span></span>
             <div style="display:flex;gap:6px;">
               <span class="site-tag variant-site-tag">📍 Límite Canónico</span>
               <span class="site-tag stop-site-tag">🛑 STOP Prematuro</span>
             </div>`,
            mutTrack
        ));
        wrap.appendChild(mutSection);
    }

    containerEl.appendChild(wrap);
}

/**
 * Asynchronously populates the SpliceAI Oracle Predictive Card
 */
async function populateSpliceAIOracle(cardElement, model) {
    const contentArea = cardElement.querySelector('#spliceaiContentArea');
    const badge = cardElement.querySelector('#spliceaiStatusBadge');
    if (!contentArea) return;

    try {
        const pred = await fetchSpliceAIPrediction(model.chromosome, model.variant.pos, model.variant.ref, model.variant.alt);
        if (!pred || !pred.available) {
            if (badge) {
                badge.className = 'badge badge-neutral';
                badge.textContent = 'Sin datos precalculados';
            }
            contentArea.innerHTML = `
                <div style="font-size: 0.88rem; color: var(--text-secondary); line-height: 1.5; padding: 4px 0;">
                    ℹ️ <em>No se encontraron puntajes precalculados de SpliceAI en GeneBE para esta coordenada genómica. Puedes continuar evaluando los escenarios de splicing mecánicos a continuación.</em>
                </div>
            `;
            return;
        }

        if (badge) {
            badge.className = `badge ${pred.recommendation?.badgeClass || 'badge-neutral'}`;
            badge.textContent = pred.recommendation?.badge || `Max DS: ${pred.maxScore.toFixed(2)}`;
        }

        const getMetricThresholdClass = (val) => {
            if (val >= 0.8) return 'threshold-very-high';
            if (val >= 0.5) return 'threshold-high';
            if (val >= 0.2) return 'threshold-suspect';
            return 'threshold-low';
        };

        const getMetricConfidenceBadge = (val) => {
            if (val >= 0.8) return '<span class="badge badge-danger" style="font-size:0.75rem;">Muy Alto (≥ 0.80)</span>';
            if (val >= 0.5) return '<span class="badge badge-warning" style="font-size:0.75rem;">Alto (≥ 0.50)</span>';
            if (val >= 0.2) return '<span class="badge badge-info" style="font-size:0.75rem;">Sospechoso (≥ 0.20)</span>';
            return '<span class="badge badge-neutral" style="font-size:0.75rem;">Bajo (&lt; 0.20)</span>';
        };

        const formatDp = (dp) => (dp > 0 ? `+${dp} pb` : `${dp} pb`);
        const isLoss = pred.recommendation?.action === 'loss';
        const isGain = pred.recommendation?.action === 'gain';

        const rows = [
            {
                name: 'Aceptor Gain (AG)',
                score: pred.ds_ag,
                dp: pred.dp_ag,
                effect: 'Ganancia de aceptor críptico (elongación / inclusión aberrante)'
            },
            {
                name: 'Aceptor Loss (AL)',
                score: pred.ds_al,
                dp: pred.dp_al,
                effect: 'Pérdida de aceptor canónico (salto de exón o retención)'
            },
            {
                name: 'Donor Gain (DG)',
                score: pred.ds_dg,
                dp: pred.dp_dg,
                effect: 'Ganancia de dador críptico (truncamiento / uso alternativo)'
            },
            {
                name: 'Donor Loss (DL)',
                score: pred.ds_dl,
                dp: pred.dp_dl,
                effect: 'Pérdida de dador canónico (salto de exón o retención)'
            }
        ];

        const tableRowsHtml = rows.map(r => {
            const isSignificant = r.score >= 0.2;
            const rowStyle = isSignificant ? 'style="background: rgba(244, 63, 94, 0.06); font-weight: 600;"' : '';
            return `
                <tr ${rowStyle}>
                    <td>
                        <strong style="color: var(--text-primary); font-size: 0.88rem;">${r.name}</strong>
                    </td>
                    <td>
                        <span class="spliceai-score-pill ${getMetricThresholdClass(r.score)}">
                            ${r.score.toFixed(3)}
                        </span>
                    </td>
                    <td style="font-family: var(--font-mono); font-size: 0.88rem; font-weight: 700; color: var(--text-primary);">
                        ${formatDp(r.dp)}
                    </td>
                    <td>
                        ${getMetricConfidenceBadge(r.score)}
                    </td>
                    <td style="font-size: 0.84rem; color: var(--text-secondary); line-height: 1.4;">
                        ${r.effect}
                    </td>
                </tr>
            `;
        }).join('');

        contentArea.innerHTML = `
            <div class="spliceai-table-wrapper">
                <table class="spliceai-table">
                    <thead>
                        <tr>
                            <th>Evento</th>
                            <th>Delta Score (DS)</th>
                            <th>Delta Posición (DP)</th>
                            <th>Nivel de Confianza</th>
                            <th>Consecuencia Prevista</th>
                        </tr>
                    </thead>
                    <tbody>
                        ${tableRowsHtml}
                    </tbody>
                </table>
            </div>

            <div class="spliceai-recommendation-box ${isLoss ? 'is-loss' : (isGain ? 'is-gain' : '')}">
                <div class="spliceai-recommendation-title">
                    ${pred.recommendation?.title || 'Evaluación de Splicing'}
                </div>
                <p class="spliceai-recommendation-text">
                    ${pred.recommendation?.text || ''}
                </p>
                ${isLoss ? `
                    <div style="margin-top: 12px; display: flex; gap: 8px; flex-wrap: wrap;">
                        <button id="btnOracleGoScenario1" class="btn-secondary" style="padding: 6px 14px; font-size: 0.82rem; cursor: pointer;">
                            Simular Escenario 1 (Salto de Exón)
                        </button>
                        <button id="btnOracleGoScenario3" class="btn-secondary" style="padding: 6px 14px; font-size: 0.82rem; cursor: pointer;">
                            Simular Escenario 3 (Retención de Intrón)
                        </button>
                    </div>
                ` : ''}
                ${isGain ? `
                    <div style="margin-top: 12px;">
                        <button id="btnOracleApplyCryptic" class="btn-secondary" style="padding: 6px 14px; font-size: 0.82rem; cursor: pointer;">
                            Configurar Escenario 2 con Δ = ${pred.recommendation?.deltaPos || 4} pb
                        </button>
                    </div>
                ` : ''}
            </div>
        `;

        // Wire quick action buttons
        const btnGoS1 = contentArea.querySelector('#btnOracleGoScenario1');
        if (btnGoS1) {
            btnGoS1.addEventListener('click', () => {
                const tabS1 = document.getElementById('tabScenario1');
                if (tabS1) tabS1.click();
                const btnRecalcS1 = document.getElementById('btnRecalcScenario1');
                if (btnRecalcS1) btnRecalcS1.click();
            });
        }

        const btnGoS3 = contentArea.querySelector('#btnOracleGoScenario3');
        if (btnGoS3) {
            btnGoS3.addEventListener('click', () => {
                const tabS3 = document.getElementById('tabScenario3');
                if (tabS3) tabS3.click();
                const btnRecalcS3 = document.getElementById('btnRecalcScenario3');
                if (btnRecalcS3) btnRecalcS3.click();
            });
        }

        const btnApplyCryptic = contentArea.querySelector('#btnOracleApplyCryptic');
        if (btnApplyCryptic) {
            btnApplyCryptic.addEventListener('click', () => {
                const tabS2 = document.getElementById('tabScenario2');
                if (tabS2) tabS2.click();
                const deltaInput = document.getElementById('crypticDeltaInput');
                const locSelect = document.getElementById('crypticLocationType');
                const dp = pred.recommendation?.deltaPos || 4;
                if (deltaInput) deltaInput.value = Math.abs(dp);
                if (locSelect) locSelect.value = dp > 0 ? 'intron' : 'exon';
                const btnRecalcS2 = document.getElementById('btnRecalcScenario2');
                if (btnRecalcS2) btnRecalcS2.click();
            });
        }

    } catch (err) {
        console.warn("Error rendering SpliceAI Oracle card:", err);
        if (contentArea) {
            contentArea.innerHTML = `
                <div style="font-size: 0.85rem; color: var(--text-muted);">
                    ℹ️ SpliceAI no disponible temporalmente en GeneBE. Puedes continuar con la simulación interactiva abajo.
                </div>
            `;
        }
    }
}

export function renderConsequenceSimulator(container, model) {
    if (!container || !model) return;
    container.innerHTML = '';

    const { exons, introns, chromosome, variant, variantLocation } = model;
    const isSplicing = variantLocation.type === 'intron' || variantLocation.isSpliceSite;
    const affectedExonNum = variantLocation.exon?.exonNum || variantLocation.intron?.donorExon || 1;
    const adjacentIntron = variantLocation.intron || introns.find(i => i.donorExon === affectedExonNum) || introns[0];

    const card = document.createElement('div');
    card.className = 'consequence-sim-card';

    if (isSplicing) {
        card.innerHTML = `
            <div class="sim-header">
                <div>
                    <h3 style="margin: 0 0 4px 0; color: var(--accent-cyan); display: flex; align-items: center; gap: 8px;">
                        <span>🔬 Simulador de Splicing</span>
                        <span class="badge badge-warning" style="font-size:0.75rem;">Variante de Splicing</span>
                    </h3>
                    <p style="margin: 0; font-size: 0.85rem; color: var(--text-secondary);">
                        Evalúa de forma interactiva y visual los 3 escenarios biológicos principales cuando se afecta un sitio de empalme.
                    </p>
                </div>
            </div>

            <div class="sim-body" style="margin-top: 16px;">
                <!-- ═══════ ORÁCULO PREDICTIVO SPLICEAI (GENEBE) ═══════ -->
                <div class="spliceai-oracle-card" id="spliceaiOracleCard">
                    <div class="spliceai-oracle-header">
                        <div class="spliceai-oracle-title">
                            <span>Predicción SpliceAI</span>
                            <span id="spliceaiStatusBadge" class="badge badge-neutral">Consultando...</span>
                        </div>
                    </div>
                    <div id="spliceaiContentArea">
                        <div style="display:flex; align-items:center; gap:10px; color:var(--text-secondary); font-size:0.9rem; padding:12px 0;">
                            <span class="calc-spinner">⏳</span> Cargando predicciones de SpliceAI...
                        </div>
                    </div>
                </div>

                <!-- Scenario Selector Tabs -->
                <div class="splicing-scenario-tabs">
                    <button class="scenario-tab-btn" id="tabScenario1" data-scenario="1">
                        1. Salto de Exón Completo (Exon Skipping)
                    </button>
                    <button class="scenario-tab-btn" id="tabScenario2" data-scenario="2">
                        2. Activación de Sitios Crípticos
                    </button>
                    <button class="scenario-tab-btn" id="tabScenario3" data-scenario="3">
                        3. Retención de Intrón Completo
                    </button>
                </div>

                <!-- Initial prompt when no scenario is selected -->
                <div id="noScenarioSelectedBox" class="nmd-decision-box" style="margin-top: 16px;">
                    <div class="nmd-title"><span>👉 Selecciona un escenario de splicing:</span></div>
                    <p class="nmd-description">Haz clic en cualquiera de los 3 botones superiores (<strong>Salto de Exón</strong>, <strong>Sitio Críptico</strong> o <strong>Retención de Intrón</strong>) para desplegar sus parámetros y simulador.</p>
                </div>

                <!-- ═══════ ESCENARIO 1: SALTO DE EXÓN COMPLETO ═══════ -->
                <div id="panelScenario1" class="splicing-scenario-panel" style="display: none;">
                    <div class="scenario-intro-box">
                        <strong>📌 Escenario 1: Salto de Exón Completo (*Exon Skipping*)</strong><br>
                        Frente a la pérdida de un sitio dador (+1/+2) o aceptor (-1/-2), el spliceosoma se saltea el exón completo y une el exón previo con el posterior.
                    </div>

                    <div style="margin-bottom: 8px; font-size: 0.85rem; font-weight: 600; color: var(--text-primary);">
                        Selecciona el exón omitido en el ARNm maduro:
                    </div>
                    <div class="exon-selector-chips" id="scenario1ExonChips">
                        ${exons.map(ex => `
                            <button class="exon-chip-btn ${ex.exonNum === affectedExonNum ? 'selected-skip' : ''}" data-exon="${ex.exonNum}">
                                Exón ${ex.exonNum} (${ex.length} pb, Fase ${ex.phase}➔${ex.endPhase})
                            </button>
                        `).join('')}
                    </div>

                    <div style="margin-top: 14px;">
                        <button id="btnRecalcScenario1" class="btn-primary">
                            ⚡ Recalcular y Comparar Pistas WT vs Salto de Exón
                        </button>
                    </div>

                    <div id="resultScenario1" class="nmd-decision-box" style="margin-top: 16px;">
                        <div class="nmd-title"><span>Evaluación Didáctica del Salto de Exón:</span></div>
                        <p class="nmd-description">Selecciona el exón y presiona calcular para comparar la versión WT frente a la unión modificada.</p>
                    </div>

                    <div id="scenario1ViewerContainer"></div>
                </div>

                <!-- ═══════ ESCENARIO 2: ACTIVACIÓN DE SITIOS CRÍPTICOS ═══════ -->
                <div id="panelScenario2" class="splicing-scenario-panel" style="display: none;">
                    <div class="scenario-intro-box">
                        <strong>📌 Escenario 2: Activación de Sitios Crípticos (*Sitios Alternativos*)</strong><br>
                        El spliceosoma utiliza una secuencia de consenso críptica cercana dentro del exón (deleción parcial) o dentro del intrón (inserción de bases intrónicas). Se visualiza en dos pasos claros: estructura primaria pre-ARNm y empalme mutado con regiones identificadas, mutación visible y traducción.
                    </div>

                    <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(260px, 1fr)); gap: 14px; margin-bottom: 14px;">
                        <div class="input-card" style="padding: 12px;">
                            <label style="font-size: 0.82rem; font-weight: 700; color: var(--accent-cyan); display: block; margin-bottom: 6px;">
                                Ubicación del Sitio Críptico:
                            </label>
                            <select id="crypticLocationType" class="input-field" style="width: 100%;">
                                <option value="intron">Dentro del Intrón (Inserción parcial de bases intrónicas)</option>
                                <option value="exon">Dentro del Exón (Deleción parcial de bases exónicas)</option>
                            </select>
                        </div>
                        <div class="input-card" style="padding: 12px;">
                            <label style="font-size: 0.82rem; font-weight: 700; color: var(--accent-cyan); display: block; margin-bottom: 6px;">
                                Distancia al Sitio Original (Δ pb):
                            </label>
                            <div style="display: flex; gap: 8px;">
                                <input type="number" id="crypticDeltaInput" class="input-field" value="4" min="1" max="500" style="flex: 1;">
                                <select id="crypticPresets" class="input-field" style="flex: 1.4;">
                                    <option value="4">+4 pb (Out-of-frame)</option>
                                    <option value="7">+7 pb (Out-of-frame)</option>
                                    <option value="12">+12 pb (In-frame, +4 aa)</option>
                                    <option value="-5">-5 pb (Out-of-frame)</option>
                                    <option value="-9">-9 pb (In-frame, -3 aa)</option>
                                </select>
                            </div>
                        </div>
                    </div>

                    <button id="btnRecalcScenario2" class="btn-primary">
                        ⚡ Recalcular y Comparar Pistas WT vs Sitio Críptico
                    </button>

                    <div id="resultScenario2" class="nmd-decision-box" style="margin-top: 16px;">
                        <div class="nmd-title"><span>Evaluación Didáctica del Sitio Críptico:</span></div>
                        <p class="nmd-description">Configura la distancia en pares de bases y presiona recalcular.</p>
                    </div>

                    <div id="scenario2ViewerContainer"></div>
                </div>

                <!-- ═══════ ESCENARIO 3: RETENCIÓN DE INTRÓN COMPLETO ═══════ -->
                <div id="panelScenario3" class="splicing-scenario-panel" style="display: none;">
                    <div class="scenario-intro-box">
                        <strong>📌 Escenario 3: Retención de Intrón (*Intron Retention*)</strong><br>
                        Fallo completo del reconocimiento de los sitios de corte; el intrón completo se retiene en el ARNm maduro y se traduce de forma continua heredando el marco del exón hasta encontrar el codón de parada prematuro (PTC) en la secuencia.
                    </div>
                    <div style="margin-bottom: 12px; font-size: 0.85rem; color: var(--text-secondary);">
                        Intrón Retenido: <strong>Intrón ${adjacentIntron ? adjacentIntron.intronNum : 1}</strong> (${adjacentIntron ? adjacentIntron.length.toLocaleString() : 0} pb).
                    </div>
                    <button id="btnRecalcScenario3" class="btn-primary">
                        ⚡ Visualizar Retención Intrónica y Búsqueda de Stop en Secuencia
                    </button>

                    <div id="resultScenario3" class="nmd-decision-box" style="margin-top: 16px;">
                        <div class="nmd-title"><span>Evaluación de Retención de Intrón:</span></div>
                        <p class="nmd-description">Presiona el botón superior para calcular y visualizar la traducción continua hacia el intrón hasta el Stop.</p>
                    </div>

                    <div id="scenario3ViewerContainer"></div>
                </div>
            </div>
        `;
    } else {
        card.innerHTML = `
            <div class="sim-header">
                <div>
                    <h3 style="margin: 0 0 4px 0; color: var(--accent-cyan); display: flex; align-items: center; gap: 8px;">
                        <span>🔬 Simulador de Marco de Lectura y Codón Stop Prematuro (PTC)</span>
                        <span class="badge ${variantLocation.type === 'exon' ? 'badge-primary' : 'badge-neutral'}">Variante Exónica</span>
                    </h3>
                    <p style="margin: 0; font-size: 0.85rem; color: var(--text-secondary);">
                        Visualiza los nuevos aminoácidos traducidos en la fase desplazada tras una inserción, deleción o sustitución.
                    </p>
                </div>
            </div>

            <div class="sim-body">
                <!-- SpliceAI card also available for exonic variants to detect exonic cryptic activation or splice disruption -->
                <div class="spliceai-oracle-card" id="spliceaiOracleCard">
                    <div class="spliceai-oracle-header">
                        <div class="spliceai-oracle-title">
                            <span>Predicción SpliceAI</span>
                            <span id="spliceaiStatusBadge" class="badge badge-neutral">Consultando...</span>
                        </div>
                    </div>
                    <div id="spliceaiContentArea">
                        <div style="display:flex; align-items:center; gap:10px; color:var(--text-secondary); font-size:0.9rem; padding:12px 0;">
                            <span class="calc-spinner">⏳</span> Cargando predicciones de SpliceAI...
                        </div>
                    </div>
                </div>

                <div style="margin-bottom: 14px;">
                    <button id="btnRunFrameshiftModule" class="btn-primary" style="padding: 12px 24px;">
                        ⚡ Calcular Desplazamiento del Marco y Visualizar Nuevos Aminoácidos
                    </button>
                </div>

                <div id="frameshiftResultBox" class="nmd-decision-box" style="margin-top: 14px;">
                    <div class="nmd-title">
                        <span>Pauta de Lectura Mutada:</span>
                    </div>
                    <div class="nmd-description" id="frameshiftResultText">
                        Presiona el botón superior para calcular los codones en la nueva fase de lectura tras la variante <strong>${variant.ref} &gt; ${variant.alt}</strong>.
                    </div>
                </div>

                <div id="frameshiftViewerContainer" style="margin-top: 18px; display: none;"></div>
            </div>
        `;
    }

    container.appendChild(card);

    // Asynchronously populate SpliceAI Oracle Card
    populateSpliceAIOracle(card, model);

    // ─── EVENT HANDLERS AND DUAL VIEWERS FOR SPLICING ───
    if (isSplicing) {
        const tabs = [
            { btn: card.querySelector('#tabScenario1'), panel: card.querySelector('#panelScenario1') },
            { btn: card.querySelector('#tabScenario2'), panel: card.querySelector('#panelScenario2') },
            { btn: card.querySelector('#tabScenario3'), panel: card.querySelector('#panelScenario3') }
        ];
        const placeholder = card.querySelector('#noScenarioSelectedBox');

        tabs.forEach(({ btn, panel }) => {
            if (btn && panel) {
                btn.addEventListener('click', () => {
                    const isAlreadyActive = btn.classList.contains('active');
                    tabs.forEach(t => {
                        t.btn?.classList.remove('active');
                        if (t.panel) t.panel.style.display = 'none';
                    });
                    if (placeholder) placeholder.style.display = 'none';

                    if (!isAlreadyActive) {
                        btn.classList.add('active');
                        panel.style.display = 'block';
                    } else if (placeholder) {
                        placeholder.style.display = 'block';
                    }
                });
            }
        });

        // Scenario 1: Exon Skipping Dual Comparison (with 5' -> 3' strand -1 support)
        let selectedExonNum = affectedExonNum;
        card.querySelectorAll('#scenario1ExonChips .exon-chip-btn').forEach(chip => {
            chip.addEventListener('click', () => {
                card.querySelectorAll('#scenario1ExonChips .exon-chip-btn').forEach(c => c.classList.remove('selected-skip'));
                chip.classList.add('selected-skip');
                selectedExonNum = parseInt(chip.dataset.exon, 10);
            });
        });

        const btnS1 = card.querySelector('#btnRecalcScenario1');
        const resS1 = card.querySelector('#resultScenario1');
        const viewerS1 = card.querySelector('#scenario1ViewerContainer');

        if (btnS1 && resS1) {
            btnS1.addEventListener('click', async () => {
                btnS1.disabled = true;
                btnS1.innerHTML = '<span class="calc-spinner">⏳ Cargando secuencias...</span>';
                if (viewerS1) viewerS1.innerHTML = '';

                let targetIdx = exons.findIndex(e => e.exonNum === selectedExonNum);
                if (targetIdx === -1) targetIdx = 0;

                const targetEx = exons[targetIdx] || exons[0];
                const prevEx = targetIdx > 0 ? exons[targetIdx - 1] : null;
                const nextEx = targetIdx < exons.length - 1 ? exons[targetIdx + 1] : null;

                const skippedLen = targetEx.codingLen || targetEx.length || 60;
                const inFrame = (skippedLen % 3 === 0);
                const shift = skippedLen % 3;

                let seqPrev = 'CAGTACCAGTTGAC';
                let seqTarget = 'TTTGAAAGTGATGAA';
                let seqTargetTail = null;
                let seqNext = 'TTAGCTGAATTGGAC';
                let seqPrevPhase = 0;
                let seqTargetPhase = targetEx.phase ?? 0;
                let seqNextPhase = nextEx ? (nextEx.phase ?? 0) : 0;

                if (prevEx) {
                    const tailInfo = getCodonAlignedExonTail(prevEx, 15);
                    const fetchedUp = await fetchExonSegment5to3(chromosome, tailInfo, model.strandNumeric);
                    if (fetchedUp) seqPrev = fetchedUp;
                    seqPrevPhase = tailInfo.phase; // 0
                }

                if (targetEx) {
                    const headInfo = getCodonAlignedExonHead(targetEx, 15);
                    const fetchTarget = await fetchExonSegment5to3(chromosome, headInfo, model.strandNumeric);
                    if (fetchTarget) seqTarget = fetchTarget;
                    seqTargetPhase = headInfo.phase;

                    if (skippedLen > 30) {
                        const tailInfo = getCodonAlignedExonTail(targetEx, 15);
                        const fetchTargetTail = await fetchExonSegment5to3(chromosome, tailInfo, model.strandNumeric);
                        if (fetchTargetTail) seqTargetTail = fetchTargetTail;
                    }
                }

                if (nextEx) {
                    const headInfo = getCodonAlignedExonHead(nextEx, 15);
                    const fetchedDown = await fetchExonSegment5to3(chromosome, headInfo, model.strandNumeric);
                    if (fetchedDown) seqNext = fetchedDown;
                    seqNextPhase = headInfo.phase;
                }

                resS1.className = `nmd-decision-box ${inFrame ? '' : 'trigger-nmd'}`;
                resS1.innerHTML = `
                    <div class="nmd-title">
                        <span>${inFrame ? '✅ Salto de Exón EN MARCO (In-Frame)' : '⚠️ Salto de Exón FUERA DE MARCO (Frameshift)'}</span>
                    </div>
                    <div class="nmd-description">
                        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; margin-bottom: 12px;">
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>1. Exón Omitido:</strong> Exón ${selectedExonNum} (${skippedLen} pb)
                            </div>
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>2. Ecuación de Marco:</strong> ${skippedLen} mod 3 = <strong>${shift}</strong>
                            </div>
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>3. Nuevo Empalme:</strong> ${prevEx ? `Exón ${prevEx.exonNum}` : 'Inicio'} ➔ ${nextEx ? `Exón ${nextEx.exonNum}` : 'Fin'}
                            </div>
                        </div>
                        ${inFrame
                            ? `&bull; <strong>Consecuencia Docente:</strong> Se eliminan exactamente <strong>${skippedLen / 3} aminoácidos</strong>. El marco de lectura aguas abajo queda <em>completamente intacto</em>.`
                            : `&bull; <strong>Consecuencia Docente:</strong> Como ${skippedLen} no es múltiplo de 3 (resto = +${shift} nt), se <em>rompe el marco de lectura</em> en la unión al ${nextEx ? `Exón ${nextEx.exonNum}` : 'siguiente exón'}. Codones aberrantes hasta <strong>🛑 STOP *</strong> prematuro.`
                        }
                    </div>
                `;

                if (viewerS1) {
                    renderComparisonSplicingViewer(viewerS1, {
                        mode: 'skipping',
                        prevExon: prevEx,
                        targetExon: targetEx,
                        nextExon: nextEx,
                        skippedLen,
                        isFrameshift: !inFrame,
                        frameshiftShift: shift,
                        seqPrev, seqPrevPhase,
                        seqTarget, seqTargetTail, seqTargetPhase,
                        seqNext, seqNextPhase
                    });
                }
                btnS1.disabled = false;
                btnS1.innerHTML = '⚡ Recalcular y Comparar Pistas WT vs Salto de Exón';
            });
        }

        // Scenario 2: Cryptic Site Dual Comparison (with strand -1 biological orientation)
        const crypticLoc = card.querySelector('#crypticLocationType');
        const crypticDelta = card.querySelector('#crypticDeltaInput');
        const crypticPresets = card.querySelector('#crypticPresets');
        const btnS2 = card.querySelector('#btnRecalcScenario2');
        const resS2 = card.querySelector('#resultScenario2');
        const viewerS2 = card.querySelector('#scenario2ViewerContainer');

        if (crypticPresets && crypticDelta && crypticLoc) {
            crypticPresets.addEventListener('change', () => {
                const val = parseInt(crypticPresets.value, 10);
                crypticDelta.value = Math.abs(val);
                crypticLoc.value = val > 0 ? 'intron' : 'exon';
            });
        }

        if (btnS2 && resS2) {
            btnS2.addEventListener('click', async () => {
                btnS2.disabled = true;
                btnS2.innerHTML = '<span class="calc-spinner">⏳ Cargando secuencias...</span>';
                if (viewerS2) viewerS2.innerHTML = '';

                const delta = parseInt(crypticDelta.value, 10) || 4;
                const isIntron = crypticLoc.value === 'intron';
                const inFrame = (delta % 3 === 0);
                const shift = delta % 3;
                const isAntisense = (model.strandNumeric === -1);

                let seqUpstream = 'CAGTACCAGTTGAC';
                let seqDownstream = 'TTAGCTGAATTGGAC';
                let crypticIntronSeq = null;
                let variantPosIndex = null;
                let seqPrevPhase = 0;
                let seqNextPhase = 0;

                const donorEx = exons.find(e => e.exonNum === affectedExonNum) || exons[0];
                const nextEx = exons.find(e => e.exonNum === donorEx.exonNum + 1) || exons[1] || null;

                if (donorEx) {
                    const tailInfo = getCodonAlignedExonTail(donorEx, 15);
                    const fetchedUp = await fetchExonSegment5to3(chromosome, tailInfo, model.strandNumeric);
                    if (fetchedUp) seqUpstream = fetchedUp;
                    seqPrevPhase = tailInfo.phase; // 0

                    if (isIntron) {
                        let realChunk = "";
                        if (!isAntisense) {
                            const intronStart = donorEx.end + 1;
                            const intronEnd = intronStart + delta - 1;
                            let raw = await fetchRegionSequence(chromosome, intronStart, intronEnd);
                            if (!raw || raw.length === 0) raw = "GTAAGTTCCAGTGAC".substring(0, delta);

                            if (variant.pos >= intronStart && variant.pos <= intronEnd) {
                                const chunkOffset = variant.pos - intronStart;
                                const arr = raw.split('');
                                arr[chunkOffset] = variant.alt;
                                raw = arr.join('');
                                variantPosIndex = seqUpstream.length + chunkOffset;
                            }
                            realChunk = raw;
                        } else {
                            // Antisense: downstream into intron is donorEx.start - 1 down to donorEx.start - delta
                            const intronStart = donorEx.start - delta;
                            const intronEnd = donorEx.start - 1;
                            let raw = await fetchRegionSequence(chromosome, intronStart, intronEnd);
                            if (!raw || raw.length === 0) raw = "GTAAGTTCCAGTGAC".substring(0, delta);

                            if (variant.pos >= intronStart && variant.pos <= intronEnd) {
                                const arr = raw.split('');
                                arr[variant.pos - intronStart] = variant.alt;
                                raw = arr.join('');
                                variantPosIndex = seqUpstream.length + (intronEnd - variant.pos);
                            }
                            realChunk = reverseComplement(raw);
                        }
                        crypticIntronSeq = realChunk;
                    } else {
                        // Exon deletion: check if variant is in deleted exon portion
                        if (!isAntisense) {
                            if (variant.pos > donorEx.end - delta && variant.pos <= donorEx.end) {
                                const delOffset = variant.pos - (donorEx.end - delta + 1);
                                variantPosIndex = Math.max(0, seqUpstream.length - delta + delOffset);
                            }
                        } else {
                            if (variant.pos >= donorEx.start && variant.pos < donorEx.start + delta) {
                                const delOffset = variant.pos - donorEx.start;
                                variantPosIndex = Math.max(0, seqUpstream.length - delta + delOffset);
                            }
                        }
                    }
                }

                if (nextEx) {
                    const headInfo = getCodonAlignedExonHead(nextEx, 15);
                    const fetchedDown = await fetchExonSegment5to3(chromosome, headInfo, model.strandNumeric);
                    if (fetchedDown) seqDownstream = fetchedDown;
                    seqNextPhase = headInfo.phase;
                }

                // Detailed variant impact explanation
                let variantNote = '';
                if (isIntron && variantPosIndex !== null) {
                    variantNote = `&bull; <strong>Variante incluida:</strong> La variante <code>Chr${chromosome}:${variant.pos} (${variant.ref}&gt;${variant.alt})</code> queda dentro de la secuencia intrónica incorporada al ARNm maduro, originando el nuevo sitio críptico.<br>`;
                } else if (!isIntron && variantPosIndex !== null) {
                    variantNote = `&bull; <strong>Variante excluida:</strong> La variante <code>Chr${chromosome}:${variant.pos} (${variant.ref}&gt;${variant.alt})</code> se ubica en el segmento exónico eliminado por el corte críptico.<br>`;
                } else {
                    variantNote = `&bull; <strong>Mecanismo:</strong> La variante <code>Chr${chromosome}:${variant.pos} (${variant.ref}&gt;${variant.alt})</code> altera el sitio dador canónico y activa el corte críptico alternativo a <strong>${delta} pb</strong>.<br>`;
                }

                resS2.className = `nmd-decision-box ${inFrame ? '' : 'trigger-nmd'}`;
                resS2.innerHTML = `
                    <div class="nmd-title">
                        <span>${inFrame ? '✅ Sitio Críptico EN MARCO (In-Frame)' : '⚠️ Sitio Críptico FUERA DE MARCO (Frameshift)'}</span>
                    </div>
                    <div class="nmd-description">
                        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; margin-bottom: 12px;">
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>1. Región Modificada:</strong> ${isIntron ? `+${delta} pb intrónicos incorporados al ARNm` : `-${delta} pb exónicos eliminados`}
                            </div>
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>2. Ecuación de Marco:</strong> ${delta} mod 3 = <strong>${shift}</strong>
                            </div>
                        </div>
                        ${variantNote}
                        ${inFrame
                            ? `&bull; <strong>Consecuencia Docente:</strong> Se ${isIntron ? `incorporan +${delta / 3}` : `pierden -${delta / 3}`} aminoácidos exactos sin alterar la fase de los codones posteriores.`
                            : `&bull; <strong>Consecuencia Docente:</strong> Ruptura de fase (+${shift} nt) → codones aberrantes hasta <strong>🛑 STOP *</strong> prematuro.`
                        }
                    </div>
                `;

                if (viewerS2) {
                    renderComparisonSplicingViewer(viewerS2, {
                        mode: 'cryptic',
                        prevExon: donorEx,
                        nextExon: nextEx,
                        crypticDelta: delta,
                        isIntronCryptic: isIntron,
                        crypticIntronSeq: crypticIntronSeq,
                        variantPosIndex: variantPosIndex,
                        isFrameshift: !inFrame,
                        frameshiftShift: shift,
                        seqPrev: seqUpstream, seqPrevPhase,
                        seqNext: seqDownstream, seqNextPhase,
                        variant: variant
                    });
                }
                btnS2.disabled = false;
                btnS2.innerHTML = '⚡ Recalcular y Comparar Pistas WT vs Sitio Críptico';
            });
        }

        // Scenario 3: Intron Retention Dual Comparison (with 5' -> 3' real STOP search)
        const btnS3 = card.querySelector('#btnRecalcScenario3');
        const resS3 = card.querySelector('#resultScenario3');
        const viewerS3 = card.querySelector('#scenario3ViewerContainer');

        if (btnS3 && resS3) {
            btnS3.addEventListener('click', async () => {
                btnS3.disabled = true;
                btnS3.innerHTML = '<span class="calc-spinner">⏳ Cargando secuencias y buscando STOP...</span>';
                if (viewerS3) viewerS3.innerHTML = '';

                const donorEx = exons.find(e => e.exonNum === affectedExonNum) || exons[0];
                const nextEx = exons.find(e => e.exonNum === donorEx.exonNum + 1) || exons[1] || null;
                const intronAfter = (adjacentIntron || introns.find(i => i.donorExon === donorEx.exonNum) || introns[0]);
                const isAntisense = (model.strandNumeric === -1);

                let exonTail = 'AGCTACCCAGTT';
                let intronStartSeq = 'GTAAGTTAGCTAATGACTTGACCA';
                let nextExHead = 'TTAGCTGAATTG';
                let seqPrevPhase = 0;
                let seqNextPhase = nextEx ? (nextEx.phase ?? 0) : 0;
                let variantPosIndex = null;

                if (donorEx) {
                    const tailInfo = getCodonAlignedExonTail(donorEx, 15);
                    const exFetch = await fetchExonSegment5to3(chromosome, tailInfo, model.strandNumeric);
                    if (exFetch) exonTail = exFetch;
                    seqPrevPhase = tailInfo.phase; // 0
                }

                if (nextEx) {
                    const headInfo = getCodonAlignedExonHead(nextEx, 15);
                    const nFetch = await fetchExonSegment5to3(chromosome, headInfo, model.strandNumeric);
                    if (nFetch) nextExHead = nFetch;
                    seqNextPhase = headInfo.phase;
                }

                let stopFound = false;
                let stopCodon = 'TGA';
                let stopDistanceInIntron = 18;
                let insertedAAs = 5;

                if (intronAfter) {
                    let intrFetch = "";
                    if (!isAntisense) {
                        const fetchEnd = Math.min(intronAfter.end, intronAfter.start + 299);
                        let raw = await fetchRegionSequence(chromosome, intronAfter.start, fetchEnd);
                        if (!raw || raw.length === 0) raw = 'GTAAGTTAGCTAATGACTTGACCA';
                        if (variant.pos >= intronAfter.start && variant.pos <= fetchEnd) {
                            const vOffset = variant.pos - intronAfter.start;
                            const arr = raw.split('');
                            arr[vOffset] = variant.alt;
                            raw = arr.join('');
                            variantPosIndex = exonTail.length + vOffset;
                        }
                        intrFetch = raw;
                    } else {
                        // Antisense: downstream from donor (intronAfter.end) downwards
                        const fetchStart = Math.max(intronAfter.start, intronAfter.end - 299);
                        let raw = await fetchRegionSequence(chromosome, fetchStart, intronAfter.end);
                        if (!raw || raw.length === 0) raw = 'GTAAGTTAGCTAATGACTTGACCA';
                        if (variant.pos >= fetchStart && variant.pos <= intronAfter.end) {
                            const arr = raw.split('');
                            arr[variant.pos - fetchStart] = variant.alt;
                            raw = arr.join('');
                            variantPosIndex = exonTail.length + (intronAfter.end - variant.pos);
                        }
                        intrFetch = reverseComplement(raw);
                    }

                    // Combined transcript: Exon Tail (Phase 0) + Intron
                    const seqCombined = exonTail + intrFetch;
                    const L_exon = exonTail.length;

                    // Scan triplet by triplet from start of Exon Tail (Phase 0)
                    for (let c = 0; c + 2 < seqCombined.length; c += 3) {
                        const codon = seqCombined.substring(c, c + 3).toUpperCase();
                        const aa = CODON_TABLE[codon];
                        // If this codon is inside the intron (or split across the junction)
                        if (c >= L_exon - (L_exon % 3)) {
                            if (aa === '*') {
                                stopFound = true;
                                const stopEndIdx = c + 3;
                                stopCodon = codon;
                                stopDistanceInIntron = stopEndIdx - L_exon;
                                intronStartSeq = seqCombined.substring(L_exon, stopEndIdx);
                                break;
                            } else {
                                insertedAAs++;
                            }
                        }
                    }

                    if (!stopFound) {
                        intronStartSeq = intrFetch.substring(0, Math.min(75, intrFetch.length));
                        stopDistanceInIntron = intronStartSeq.length;
                    }
                }

                resS3.className = 'nmd-decision-box trigger-nmd';
                resS3.innerHTML = `
                    <div class="nmd-title">
                        <span>⚠️ Truncación Prematura por Retención de Intrón ${intronAfter ? intronAfter.intronNum : 1}:</span>
                    </div>
                    <div class="nmd-description">
                        <strong>📌 Foco Molecular:</strong> El fallo del sitio dador retiene la secuencia del intrón a continuación del Exón ${donorEx.exonNum}.<br>
                        &bull; <strong>Variante / Sitio Afectado:</strong> <code>Chr${chromosome}:${variant.pos.toLocaleString()} (${variant.ref} &gt; ${variant.alt})</code>.<br>
                        &bull; <strong>Codón Stop Prematuro (PTC):</strong> Se identificó el codón de parada <code>${stopCodon}</code> a <strong>+${stopDistanceInIntron} pb</strong> dentro del intrón (${insertedAAs} aminoácidos nuevos traducidos antes de STOP).<br>
                        &bull; <strong>Destino:</strong> Síntesis de proteína severamente truncada y activación de degradación por NMD (*Nonsense-Mediated Decay*).
                    </div>
                `;

                if (viewerS3) {
                    renderComparisonSplicingViewer(viewerS3, {
                        mode: 'retention',
                        prevExon: donorEx,
                        nextExon: nextEx,
                        seqPrev: exonTail, seqPrevPhase,
                        seqNext: nextExHead, seqNextPhase,
                        retainedIntronSeq: intronStartSeq,
                        variantPosIndex: variantPosIndex,
                        variant: variant
                    });
                }
                btnS3.disabled = false;
                btnS3.innerHTML = '⚡ Visualizar Retención Intrónica y Búsqueda de Stop en Secuencia';
            });
        }
    } else {
        // Scenario: Frameshift Live Recalculation (with strand -1 biological 5' -> 3' support)
        const btnRunFs = card.querySelector('#btnRunFrameshiftModule');
        const fsBox = card.querySelector('#frameshiftResultBox');
        const fsViewer = card.querySelector('#frameshiftViewerContainer');

        if (btnRunFs && fsBox) {
            btnRunFs.addEventListener('click', async () => {
                btnRunFs.disabled = true;
                btnRunFs.innerHTML = '<span class="calc-spinner">⏳ Cargando secuencias y calculando marco...</span>';
                if (fsViewer) fsViewer.innerHTML = '';

                const targetExon = variantLocation.exon || exons.find(e => variant.pos >= e.start && variant.pos <= e.end) || exons[0];
                const exonEntryPhase = targetExon.phase ?? 0;
                const isAntisense = (model.strandNumeric === -1);

                let fullSeq = "";
                let relVarIdx = 0;
                let startPhase = 0;

                if (!isAntisense) {
                    const varOffsetInExon = Math.max(0, variant.pos - targetExon.start);
                    const codonLead = (varOffsetInExon + exonEntryPhase) % 3;
                    const codonStartPos = variant.pos - codonLead;

                    let viewStart = Math.max(targetExon.start, codonStartPos - 15);
                    const rem = (viewStart - targetExon.start + exonEntryPhase) % 3;
                    if (rem !== 0) viewStart -= rem;
                    if (viewStart < targetExon.start) viewStart = targetExon.start;
                    startPhase = (viewStart - targetExon.start + exonEntryPhase) % 3; // 0

                    const fetchEnd = Math.min(model.end, Math.max(targetExon.end, variant.pos + 300));
                    let raw = await fetchRegionSequence(chromosome, viewStart, fetchEnd);
                    if (!raw || raw.length === 0) raw = "N".repeat(fetchEnd - viewStart + 1);

                    fullSeq = raw;
                    relVarIdx = variant.pos - viewStart;
                } else {
                    // Antisense: 5' entrance is targetExon.end, moving down to targetExon.start
                    const varOffsetInExon = Math.max(0, targetExon.end - variant.pos);
                    const codonLead = (varOffsetInExon + exonEntryPhase) % 3;
                    const codonStartPos = variant.pos + codonLead; // higher genomic coord

                    let viewStartCoord = Math.min(targetExon.end, codonStartPos + 15);
                    const rem = (targetExon.end - viewStartCoord + exonEntryPhase) % 3;
                    if (rem !== 0) viewStartCoord += rem;
                    if (viewStartCoord > targetExon.end) viewStartCoord = targetExon.end;
                    startPhase = (targetExon.end - viewStartCoord + exonEntryPhase) % 3; // 0

                    const fetchEndCoord = Math.max(model.start, Math.min(targetExon.start, variant.pos - 300));
                    let rawGenomic = await fetchRegionSequence(chromosome, fetchEndCoord, viewStartCoord);
                    if (!rawGenomic || rawGenomic.length === 0) rawGenomic = "N".repeat(viewStartCoord - fetchEndCoord + 1);

                    fullSeq = reverseComplement(rawGenomic);
                    relVarIdx = viewStartCoord - variant.pos;
                }

                const effectiveRef = isAntisense ? reverseComplement(variant.ref || 'N') : (variant.ref || 'N');
                const effectiveAlt = isAntisense ? reverseComplement(variant.alt || 'N') : (variant.alt || 'N');
                const refLen = Math.max(1, effectiveRef.length);
                const deltaNt = effectiveAlt.length - effectiveRef.length;
                const shift = (deltaNt % 3 + 3) % 3;
                const isFrameshift = (shift !== 0);

                // Construct WT sequence and Mutant sequence in 5' -> 3' mRNA
                const wtSeq = fullSeq;
                const mutSeq = fullSeq.substring(0, relVarIdx) + effectiveAlt + fullSeq.substring(relVarIdx + refLen);

                // Scan mutant sequence for the first in-frame STOP codon
                let stopFound = false;
                let stopCodon = '';
                let stopEndIdx = null;
                let stopDistNt = 0;
                let newAAs = 0;

                const scanStart = startPhase === 0 ? 0 : (3 - startPhase);
                for (let c = scanStart; c + 2 < mutSeq.length; c += 3) {
                    if (c >= relVarIdx) {
                        const codon = mutSeq.substring(c, c + 3).toUpperCase();
                        const aa = CODON_TABLE[codon];
                        if (aa === '*') {
                            stopFound = true;
                            stopCodon = codon;
                            stopEndIdx = c + 3;
                            stopDistNt = stopEndIdx - relVarIdx;
                            break;
                        } else {
                            newAAs++;
                        }
                    }
                }

                // Window slice for display: up to the stop codon or up to ~75 pb
                const displayLen = stopEndIdx ? Math.min(stopEndIdx + 9, mutSeq.length) : Math.min(mutSeq.length, relVarIdx + 75);
                const mutDisplaySeq = mutSeq.substring(0, displayLen);
                const wtDisplaySeq = wtSeq.substring(0, Math.min(wtSeq.length, displayLen));

                fsBox.className = `nmd-decision-box ${isFrameshift ? 'trigger-nmd' : ''}`;
                fsBox.innerHTML = `
                    <div class="nmd-title">
                        <span>${isFrameshift ? '⚠️ Ruptura del Marco de Lectura (Frameshift Activo)' : '✅ Conservación del Marco de Lectura (In-Frame)'}</span>
                    </div>
                    <div class="nmd-description">
                        <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 10px; margin-bottom: 12px;">
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>1. Exón Afectado:</strong> Exón ${targetExon.exonNum} (Fase entrada: ${targetExon.phase} ➔ salida: ${targetExon.endPhase})
                            </div>
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>2. Locus y Cambio:</strong> Chr${chromosome}:${variant.pos.toLocaleString()} (${variant.ref} &gt; ${variant.alt}) [Δ ${deltaNt >= 0 ? `+${deltaNt}` : `${deltaNt}`} pb]
                            </div>
                            <div style="background: rgba(255,255,255,0.05); padding: 8px 12px; border-radius: 4px;">
                                <strong>3. Ecuación de Marco:</strong> ${deltaNt} mod 3 = <strong>${shift} nt</strong> (${isFrameshift ? `Desplazamiento +${shift}` : 'In-Frame'})
                            </div>
                        </div>
                        ${isFrameshift 
                            ? `&bull; <strong>Consecuencia Molecular:</strong> Ruptura de la pauta de lectura respetando el marco nativo del Exón ${targetExon.exonNum}. ${stopFound ? `Se identificó el codón de parada prematuro (PTC) <strong>${stopCodon}</strong> a <strong>+${stopDistNt} pb</strong> (${newAAs} aminoácidos nuevos traducidos antes de la interrupción 🛑).` : 'Codones aberrantes downstream en la fase desplazada.'}<br>&bull; <strong>Destino Proteico:</strong> Síntesis de proteína truncada y degradación del transcripto mediada por <em>Nonsense-Mediated Decay (NMD)</em>.`
                            : `&bull; <strong>Consecuencia Molecular:</strong> La alteración conserva la pauta de lectura de los codones posteriores (Fase intacta). ${deltaNt === 0 ? 'Sustitución de aminoácido (missense/sinónima).' : `Inserción/Deleción de ${Math.abs(deltaNt / 3)} aminoácido(s) in-frame.`}`
                        }
                    </div>
                `;

                if (fsViewer) {
                    fsViewer.style.display = 'block';
                    fsViewer.innerHTML = '';

                    const wrap = document.createElement('div');
                    wrap.className = 'splicing-comparison-container';

                    // Track 1: WT (Strictly using Step 3 exon frame)
                    const wtTrack = renderSynchronizedGridTrack([
                        { seq: wtDisplaySeq, label: `Exón ${targetExon.exonNum} WT`, type: 'exon', startPhase }
                    ], {
                        showAminoAcids: true,
                        junctionIndexes: [relVarIdx]
                    });

                    const wtBox = document.createElement('div');
                    wtBox.className = 'track-named-box';
                    wtBox.innerHTML = `
                        <div class="track-box-header">
                            <span>🧬 Secuencia Salvaje (WT) — Exón <strong>${targetExon.exonNum}</strong></span>
                            <span class="phase-tag">Fase entrada: ${targetExon.phase} ➔ Fase salida: ${targetExon.endPhase}</span>
                        </div>
                    `;
                    wtBox.appendChild(wtTrack);
                    wrap.appendChild(wtBox);

                    // Track 2: Mutated
                    const mutTrack = renderSynchronizedGridTrack([
                        { seq: mutDisplaySeq, label: `Exón ${targetExon.exonNum} Mutado`, type: 'exon', startPhase }
                    ], {
                        showAminoAcids: true,
                        junctionIndexes: [relVarIdx],
                        isFrameshift,
                        frameshiftStartIndex: relVarIdx,
                        variantPosIndex: relVarIdx,
                        variantRef: effectiveRef,
                        variantAlt: effectiveAlt,
                        stopAtCodonStop: isFrameshift
                    });

                    const mutBox = document.createElement('div');
                    mutBox.className = 'track-named-box is-cryptic-box';
                    mutBox.innerHTML = `
                        <div class="track-box-header">
                            <span>⚡ Secuencia Reconstruida — Exón <strong>${targetExon.exonNum}</strong> (${variant.ref} &gt; ${variant.alt})</span>
                            <span class="badge ${isFrameshift ? 'badge-danger' : 'badge-neutral'}">${isFrameshift ? '⚠️ Frameshift Activo' : '✅ In-Frame'}</span>
                        </div>
                    `;
                    mutBox.appendChild(mutTrack);
                    wrap.appendChild(mutBox);

                    fsViewer.appendChild(wrap);
                }

                btnRunFs.disabled = false;
                btnRunFs.innerHTML = '⚡ Calcular Desplazamiento del Marco y Visualizar Nuevos Aminoácidos';
            });
        }
    }
}
