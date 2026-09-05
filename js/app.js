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
const variantTypeSelect = $('variantTypeSelect');
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
function init() {
    setupEventListeners();
    
    // Auto-run analysis with default example
    runAnalysis("5:169670598 G>A", "NM_004946.3");
}

function setupEventListeners() {
    if (calculateBtn) {
        calculateBtn.addEventListener('click', handleCalculateClick);
    }

    if (variantInput) {
        variantInput.addEventListener('keydown', (e) => {
            if (e.key === 'Enter') handleCalculateClick();
        });
        variantInput.addEventListener('input', hideValidationError);
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

    // Modal Acerca de (Pop-up Overlay)
    window.openAboutModal = (e) => {
        if (e && e.preventDefault) e.preventDefault();
        const modal = document.getElementById('aboutModal');
        if (modal) {
            if (typeof modal.showModal === 'function') {
                if (!modal.open) modal.showModal();
            } else {
                modal.setAttribute('open', '');
                modal.style.display = 'flex';
            }
        }
    };

    window.closeAboutModal = (e) => {
        if (e && e.preventDefault) e.preventDefault();
        const modal = document.getElementById('aboutModal');
        if (modal) {
            if (typeof modal.close === 'function') {
                if (modal.open) modal.close();
            } else {
                modal.removeAttribute('open');
                modal.style.display = 'none';
            }
        }
    };

    if (aboutBtn) {
        aboutBtn.addEventListener('click', window.openAboutModal);
    }

    if (closeAboutModalBtn) {
        closeAboutModalBtn.addEventListener('click', window.closeAboutModal);
    }

    if (aboutModal) {
        aboutModal.addEventListener('click', (e) => {
            if (e.target === aboutModal) {
                window.closeAboutModal(e);
            }
        });
    }

    window.addEventListener('keydown', (e) => {
        const modal = document.getElementById('aboutModal');
        if (e.key === 'Escape' && modal && modal.open) {
            window.closeAboutModal(e);
        }
    });

    if (copyAliasMpBtn) {
        copyAliasMpBtn.addEventListener('click', handleCopyAlias);
    }
}

function handleCalculateClick() {
    const vInput = variantInput?.value.trim();
    const txInput = transcriptInput?.value.trim();
    const typeInput = variantTypeSelect?.value;

    const validation = validateGenomicInput(vInput);
    if (!validation.valid) {
        showValidationError(validation.error);
        return;
    }
    hideValidationError();
    runAnalysis(vInput, txInput, typeInput);
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
    const alias = "lorenzo.bioinf.egc";
    try {
        await navigator.clipboard.writeText(alias);
    } catch (err) {
        prompt("Copia el Alias de Mercado Pago:", alias);
    }
    if (copyAliasSuccessMsg) {
        copyAliasSuccessMsg.style.display = 'block';
        setTimeout(() => {
            copyAliasSuccessMsg.style.display = 'none';
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
    const { variant, variantLocation, geneName, chromosome } = model;

    let category = "Variante Genómica";
    let subText = "Deducción automática";
    let badgeClass = "badge-warning";
    let desc = "";

    if (variantLocation.type === 'intron') {
        if (variantLocation.isCanonicalSplice) {
            category = "Splicing Canónico";
            subText = `Sitio crítico (${variantLocation.offset})`;
            badgeClass = "badge-danger";
            desc = `Variante en unión de empalme canónica dinucleótido GT/AG (${variantLocation.name}, offset ${variantLocation.offset}).`;
        } else if (variantLocation.isSpliceSite) {
            category = "Splicing Intrónico";
            subText = `Región yuxta-exónica (${variantLocation.offset})`;
            badgeClass = "badge-warning";
            desc = `Variante intrónica cercana (${variantLocation.offset} pb de la unión exón-intrón).`;
        } else {
            category = "Intrónica";
            subText = `Profunda (${variantLocation.offset} pb)`;
            badgeClass = "badge-neutral";
            desc = `Variante en región intrónica profunda.`;
        }
    } else if (variantLocation.type === 'exon') {
        category = "Exónica";
        subText = `${variantLocation.name}`;
        badgeClass = "badge-danger";
        desc = `Variante dentro de la secuencia codificante del ${variantLocation.name}.`;
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
    if (consequenceDescText) {
        consequenceDescText.textContent = desc;
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
        calculateBtn.innerHTML = isLoading 
            ? `<span>⏳ Cargando datos desde Ensembl...</span>` 
            : `<span>⚡ Cargar y Visualizar</span>`;
    }
}

// Start
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
} else {
    init();
}
