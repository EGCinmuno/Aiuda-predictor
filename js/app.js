/**
 * NonTool App Orchestrator
 * Integrates Ensembl REST API downloads mirroring download.py,
 * manages loading states, and coordinates Step 3 (Macro Exon/Intron Map),
 * Step 4 (Micro Base Viewer), Step 5 (Classification), and Step 6 (Consequence Simulator).
 */

import { parseGenomicInput, downloadGenomicStructure, validateGenomicInput } from './ensemblService.js';
import { renderGenomicExonMap, renderFranklinBaseViewer, renderConsequenceSimulator } from './visualRenderer.js';

// Global App State
const state = {
    currentModel: null,
    theme: 'dark',
    isLoading: false
};

// DOM References
const $ = id => document.getElementById(id);

const variantInput = $('variantInput');
const transcriptInput = $('transcriptInput');
const calculateBtn = $('calculateBtn');
const presetPills = $('presetPills');
const themeToggleBtn = $('themeToggleBtn');
const aboutBtn = $('aboutBtn');
const aboutModal = $('aboutModal');
const closeAboutModalBtn = $('closeAboutModalBtn');
const copyAliasMpBtn = $('copyAliasMpBtn');
const copyAliasSuccessMsg = $('copyAliasSuccessMsg');
const validationErrorBox = $('validationErrorBox');
const validationErrorText = $('validationErrorText');

// Display Containers
const dashboardGrid = $('dashboardGrid');
const geneTitleDisplay = $('geneTitleDisplay');
const geneSubDetailsDisplay = $('geneSubDetailsDisplay');
const variantBadgeTop = $('variantBadgeTop');
const exonMapContainer = $('exonMapContainer');
const genooxViewerContainer = $('genooxViewerContainer');
const consequenceSimContainer = $('consequenceSimContainer');
const downloadStatusBox = $('downloadStatusBox');
const downloadStatusText = $('downloadStatusText');

// Consequence / Metric Displays
const variantTypeVal = $('variantTypeVal');
const variantTypeSub = $('variantTypeSub');
const genomicCoordVal = $('genomicCoordVal');
const genomicCoordSub = $('genomicCoordSub');
const locationVal = $('locationVal');
const locationSub = $('locationSub');
const splicingImpactVal = $('splicingImpactVal');
const splicingImpactSub = $('splicingImpactSub');
const consequenceCategoryBadge = $('consequenceCategoryBadge');
const consequenceDescText = $('consequenceDescText');

/**
 * Initialization
 */
// ─── i18n ────────────────────────────────────────────────────────────────────
const TRANSLATIONS = {
    es: {
        'lang-intro':       'Plataforma pedagógica para el análisis molecular de variantes de splicing. Modela la estructura exón-intrón en coordenadas GRCh38, visualiza uniones canónicas a resolución nucleotídica y simula el impacto en el marco de lectura del ARNm.',
        'lang-step1-label': 'Variante Genómica',
        'lang-step1-fmt':   '(Chr:Pos Ref>Alt)',
        'lang-step1-sub':   'Coordenada genómica (Chr:Pos Ref>Alt)',
        'lang-step2-label': 'Transcripto',
        'lang-step2-sub':   'Identificador opcional (RefSeq / Ensembl)',
        'lang-calc-btn':    '⚡ Cargar y Visualizar',
        'lang-calc-loading':'⏳ Cargando datos desde Ensembl...',
        'lang-examples-btn':'💡 Ejemplos',
        'lang-step3-title': 'Mapa de Exones e Intrones',
        'lang-step3-sub':   'Distribución genómica del transcripto.',
        'lang-step4-title': 'Visor a Nivel de Bases y Splicing',
        'lang-step4-sub':   'Regla genómica, pista de aminoácidos, uniones canónicas y marcador de variante.',
        'lang-step5-title': 'Clasificación del Tipo de Variante',
        'lang-footer-warn': '⚠️ Aviso Importante: Este asistente es una herramienta educativa de apoyo y no reemplaza el criterio clínico. No cuenta con validación clínica para diagnóstico. Toda clasificación debe ser validada por un profesional habilitado.',
        'lang-about-btn':   'ℹ️ Acerca de',
        'lang-theme-btn':   '🌓 Tema',
        'lang-lang-btn':    'EN',
    },
    en: {
        'lang-intro':       'Educational platform for molecular analysis of splicing variants. Models exon-intron structure in GRCh38 coordinates, visualizes canonical junctions at nucleotide resolution, and simulates the impact on the mRNA reading frame.',
        'lang-step1-label': 'Genomic Variant',
        'lang-step1-fmt':   '(Chr:Pos Ref>Alt)',
        'lang-step1-sub':   'Genomic coordinate (Chr:Pos Ref>Alt)',
        'lang-step2-label': 'Transcript',
        'lang-step2-sub':   'Optional identifier (RefSeq / Ensembl)',
        'lang-calc-btn':    '⚡ Load & Visualize',
        'lang-calc-loading':'⏳ Loading data from Ensembl...',
        'lang-examples-btn':'💡 Examples',
        'lang-step3-title': 'Exon & Intron Map',
        'lang-step3-sub':   'Genomic distribution of the transcript.',
        'lang-step4-title': 'Base-Level & Splicing Viewer',
        'lang-step4-sub':   'Genomic ruler, amino acid track, canonical junctions and variant marker.',
        'lang-step5-title': 'Variant Type Classification',
        'lang-footer-warn': '⚠️ Important Notice: This tool is for educational purposes only and does not replace clinical judgment. It has not been clinically validated for diagnosis. All classifications must be validated by a qualified professional.',
        'lang-about-btn':   'ℹ️ About',
        'lang-theme-btn':   '🌓 Theme',
        'lang-lang-btn':    'ES',
    }
};

let currentLang = localStorage.getItem('spliceecgenio-lang') || 'es';

function applyLang(lang) {
    currentLang = lang;
    localStorage.setItem('spliceecgenio-lang', lang);
    const t = TRANSLATIONS[lang];
    document.querySelectorAll('[data-i18n]').forEach(el => {
        const key = el.getAttribute('data-i18n');
        if (t[key] !== undefined) el.textContent = t[key];
    });
    // Update calculate button (preserves loading state awareness)
    if (calculateBtn && !state.isLoading) {
        calculateBtn.innerHTML = `<span data-i18n="lang-calc-btn">${t['lang-calc-btn']}</span>`;
    }
    // Update lang toggle button label
    const langBtn = $('langToggleBtn');
    if (langBtn) langBtn.innerHTML = `🌐 ${t['lang-lang-btn']}`;
}
// ─────────────────────────────────────────────────────────────────────────────

function init() {
    setupEventListeners();
    applyLang(currentLang);
    // No auto-run: user must click the button
}

function setupEventListeners() {
    if (calculateBtn) {
        calculateBtn.addEventListener('click', handleCalculateClick);
    }

    if (variantInput) {
        variantInput.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') handleCalculateClick();
        });
        variantInput.addEventListener('input', () => {
            hideValidationError();
            // Clear transcript field when user changes the variant input
            if (transcriptInput) transcriptInput.value = '';
        });
    }

    // Preset pills: populate inputs WITHOUT auto-running analysis
    if (presetPills) {
        presetPills.addEventListener('click', (e) => {
            const pill = e.target.closest('.preset-pill');
            if (!pill) return;
            presetPills.querySelectorAll('.preset-pill').forEach(p => p.classList.remove('active'));
            pill.classList.add('active');

            const varCoord = pill.dataset.variant;
            const txId = pill.dataset.transcript || "";
            if (variantInput) variantInput.value = varCoord;
            if (transcriptInput) transcriptInput.value = txId;

            hideValidationError();
            // User requested: Do NOT run analysis automatically on preset click
        });
    }

    // Collapsible Examples Toggle
    const togglePresetsBtn = $('togglePresetsBtn');
    const presetsArrow = $('presetsArrow');
    if (togglePresetsBtn && presetPills) {
        togglePresetsBtn.addEventListener('click', () => {
            const isHidden = presetPills.style.display === 'none';
            presetPills.style.display = isHidden ? 'flex' : 'none';
            if (presetsArrow) presetsArrow.textContent = isHidden ? '▴' : '▾';
        });
    }

    if (themeToggleBtn) {
        themeToggleBtn.addEventListener('click', toggleTheme);
    }

    const langToggleBtn = $('langToggleBtn');
    if (langToggleBtn) {
        langToggleBtn.addEventListener('click', () => {
            applyLang(currentLang === 'es' ? 'en' : 'es');
        });
    }

    // Modal Acerca de (Pop-up Overlay)
    window.openAboutModal = (e) => {
        if (e && e.preventDefault) e.preventDefault();
        const modal = document.getElementById('aboutModal');
        if (modal) {
            modal.classList.add('active');
            modal.style.display = 'flex';
        }
    };

    window.closeAboutModal = (e) => {
        if (e && e.preventDefault) e.preventDefault();
        const modal = document.getElementById('aboutModal');
        if (modal) {
            modal.classList.remove('active');
            modal.style.display = 'none';
        }
    };

    if (aboutBtn) {
        aboutBtn.addEventListener('click', window.openAboutModal);
    }

    const closeBtn = document.getElementById('closeAboutModal') || document.getElementById('closeAboutModalBtn');
    if (closeBtn) {
        closeBtn.addEventListener('click', window.closeAboutModal);
    }

    const modalOverlay = document.getElementById('aboutModal');
    if (modalOverlay) {
        modalOverlay.addEventListener('click', (e) => {
            if (e.target === modalOverlay) {
                window.closeAboutModal(e);
            }
        });
    }

    window.addEventListener('keydown', (e) => {
        const modal = document.getElementById('aboutModal');
        if (e.key === 'Escape' && modal && (modal.classList.contains('active') || modal.style.display === 'flex')) {
            window.closeAboutModal(e);
        }
    });

    const copyBtn = document.getElementById('copyAliasBtn') || document.getElementById('copyAliasMpBtn');
    if (copyBtn) {
        copyBtn.addEventListener('click', handleCopyAlias);
    }
}

function handleCalculateClick() {
    const vInput = variantInput?.value.trim();
    const txInput = transcriptInput?.value.trim();

    const validation = validateGenomicInput(vInput);
    if (!validation.valid) {
        showValidationError(validation.error);
        return;
    }
    hideValidationError();
    runAnalysis(vInput, txInput);
}

function showValidationError(msg) {
    if (validationErrorBox && validationErrorText) {
        validationErrorText.textContent = msg;
        validationErrorBox.style.display = 'flex';
    } else {
        alert(msg);
    }
}

function hideValidationError() {
    if (validationErrorBox) {
        validationErrorBox.style.display = 'none';
    }
}

async function handleCopyAlias() {
    const alias = "lorenzo.erra.mp";
    try {
        await navigator.clipboard.writeText(alias);
    } catch (err) {
        prompt("Copia el Alias de Mercado Pago:", alias);
    }
    const successMsg = document.getElementById('copyAliasSuccessMsg');
    if (successMsg) {
        successMsg.style.display = 'block';
        setTimeout(() => {
            successMsg.style.display = 'none';
        }, 3500);
    }
}

function toggleTheme() {
    state.theme = state.theme === 'dark' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', state.theme);
    if (state.currentModel) {
        renderAll(state.currentModel);
    }
}

/**
 * Main Analysis Pipeline matching download.py
 */
async function runAnalysis(varString, txString = null, typeString = null) {
    if (state.isLoading) return;
    setLoading(true);

    // Hide steps 3+ while downloading/calculating
    if (dashboardGrid) {
        dashboardGrid.style.display = 'none';
    }

    try {
        updateStatus(`Analizando formato de variante "${varString}"...`);
        const parsed = parseGenomicInput(varString, txString, typeString);

        updateStatus(`Conectando con Ensembl REST API para Chr${parsed.chrom}:${parsed.pos}...`);
        const model = await downloadGenomicStructure(parsed, (msg) => updateStatus(msg));

        state.currentModel = model;
        updateStatus(`✅ ¡Descarga y estructura completada para ${model.geneName}!`, true);

        // Show dashboard steps
        if (dashboardGrid) {
            dashboardGrid.style.display = 'flex';
        }

        await renderAll(model);
    } catch (err) {
        console.error("Error en análisis:", err);
        alert(`❌ ${err.message}`);
        updateStatus(`⚠️ Error: ${err.message}`, false);
    } finally {
        setLoading(false);
    }
}

/**
 * Renders all dashboard panels
 */
async function renderAll(model) {
    const { geneName, transcriptId, chromosome, strand, start, end, exons, introns, variant, variantLocation } = model;

    // 1. Header Info
    if (geneTitleDisplay) {
        geneTitleDisplay.innerHTML = `<span class="step-badge">Paso 3</span> ${geneName} (${transcriptId}) — Mapa de Exones e Intrones`;
    }
    if (geneSubDetailsDisplay) {
        geneSubDetailsDisplay.textContent = `Chr ${chromosome} (${strand}) | Coordenadas Genómicas: ${start.toLocaleString()}..${end.toLocaleString()} (${(end - start + 1).toLocaleString()} pb) | ${exons.length} Exones | ${introns.length} Intrones`;
    }
    if (variantBadgeTop) {
        variantBadgeTop.innerHTML = `<span class="badge ${variantLocation.isSpliceSite ? 'badge-warning' : 'badge-danger'}">Chr${chromosome}:${variant.pos} ${variant.ref}>${variant.alt}</span>`;
    }

    // 2. Paso 3: Mapa General de Exones e Intrones
    if (exonMapContainer) {
        renderGenomicExonMap(exonMapContainer, model, async (targetPos) => {
            if (genooxViewerContainer) {
                await renderFranklinBaseViewer(genooxViewerContainer, model, targetPos);
                genooxViewerContainer.scrollIntoView({ behavior: 'smooth', block: 'center' });
            }
        });
    }

    // 3. Paso 4: Visor a Nivel de Bases (WT default)
    if (genooxViewerContainer) {
        await renderFranklinBaseViewer(genooxViewerContainer, model, variant.pos);
    }

    // 4. Paso 5: Clasificación y Métricas
    updateClassificationMetrics(model);

    // 5. Paso 6: Simulador de Consecuencia y Recálculo de Marco
    if (consequenceSimContainer) {
        renderConsequenceSimulator(consequenceSimContainer, model);
    }
}

/**
 * Updates consequence cards and badges
 */
function updateClassificationMetrics(model) {
    const { variant, variantLocation, geneName, chromosome, transcriptId } = model;

    let category = "Variante Genómica";
    let subText = "Deducción automática";
    let badgeClass = "badge-warning";

    if (variantLocation.type === 'intron') {
        if (variantLocation.isCanonicalSplice) {
            category = "Splicing Canónico";
            subText = `Sitio crítico (${variantLocation.offset})`;
            badgeClass = "badge-danger";
        } else if (variantLocation.isSpliceSite) {
            category = "Splicing Intrónico";
            subText = `Región yuxta-exónica (${variantLocation.offset})`;
            badgeClass = "badge-warning";
        } else {
            category = "Intrónica";
            subText = `Profunda (${variantLocation.offset} pb)`;
            badgeClass = "badge-neutral";
        }
    } else if (variantLocation.type === 'exon') {
        category = "Exónica";
        subText = `${variantLocation.name}`;
        badgeClass = "badge-danger";
    }

    if (variantTypeVal) variantTypeVal.textContent = category;
    if (variantTypeSub) variantTypeSub.textContent = subText;

    if (genomicCoordVal) genomicCoordVal.textContent = `Chr${chromosome}:${variant.pos.toLocaleString()}`;
    if (genomicCoordSub) genomicCoordSub.textContent = `Cambio ${variant.ref} > ${variant.alt}`;

    if (locationVal) locationVal.textContent = variantLocation.name;
    if (locationSub) locationSub.textContent = `En gen ${geneName}`;

    if (splicingImpactVal) {
        splicingImpactVal.textContent = variantLocation.isSpliceSite ? 'Afecta Splice Site' : 'Conservado';
        splicingImpactVal.style.color = variantLocation.isSpliceSite ? 'var(--accent-rose)' : 'var(--accent-emerald)';
    }
    if (splicingImpactSub) {
        splicingImpactSub.textContent = variantLocation.offset ? `Offset: ${variantLocation.offset}` : 'En exón';
    }

    if (consequenceCategoryBadge) {
        consequenceCategoryBadge.innerHTML = `<span class="badge ${badgeClass}">${category}</span>`;
    }

    const locText = variantLocation.offset ? `${variantLocation.name} (${variantLocation.offset})` : variantLocation.name;
    const txLabel = transcriptId ? ` (${transcriptId})` : '';
    const desc = `Chr${chromosome}:${variant.pos.toLocaleString()} ${variant.ref} > ${variant.alt}${txLabel} es una Variante <strong>${category}</strong> ${subText}, se localiza en el <strong>${locText}</strong> del gen <strong>${geneName}</strong>.`;

    if (consequenceDescText) {
        consequenceDescText.innerHTML = desc;
    }
}

function updateStatus(text, isSuccess = null) {
    if (!downloadStatusBox || !downloadStatusText) return;
    downloadStatusBox.style.display = 'flex';
    downloadStatusText.textContent = text;
    if (isSuccess === true) {
        downloadStatusBox.style.borderLeftColor = 'var(--accent-emerald)';
    } else if (isSuccess === false) {
        downloadStatusBox.style.borderLeftColor = 'var(--accent-rose)';
    } else {
        downloadStatusBox.style.borderLeftColor = 'var(--accent-cyan)';
    }
}

function setLoading(isLoading) {
    state.isLoading = isLoading;
    if (calculateBtn) {
        calculateBtn.disabled = isLoading;
        const t = TRANSLATIONS[currentLang];
        calculateBtn.innerHTML = isLoading
            ? `<span>${t['lang-calc-loading']}</span>`
            : `<span data-i18n="lang-calc-btn">${t['lang-calc-btn']}</span>`;
    }
}

// Start
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
} else {
    init();
}
