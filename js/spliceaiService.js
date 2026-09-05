/**
 * GeneBE SpliceAI Predictive Oracle Service
 * 
 * Fetches precomputed SpliceAI scores from GeneBE REST API:
 * 1. Primary: https://api.genebe.net/shield/v1/variant/hg38/{chr}-{pos}-{ref}-{alt}
 * 2. Fallback: https://api.genebe.net/cloud/api-public/v1/variant?chr={chr}&pos={pos}&ref={ref}&alt={alt}&genome=hg38
 * 
 * Extracts the 4 Delta Scores (DS) and 4 Delta Positions (DP):
 * - ds_ag: Acceptor Gain
 * - ds_al: Acceptor Loss
 * - ds_dg: Donor Gain
 * - ds_dl: Donor Loss
 * - dp_ag, dp_al, dp_dg, dp_dl: Respective positions in bp relative to variant
 */

const GENEBE_API_URL = "https://api.genebe.net";

/**
 * Queries SpliceAI predictions for a variant in GRCh38
 * @param {string|number} chrom Chromosome (e.g. "5" or 5)
 * @param {number} pos 1-based genomic position
 * @param {string} ref Reference allele
 * @param {string} alt Alternative allele
 * @param {function} onProgress Progress callback
 * @returns {Promise<Object>} SpliceAI prediction model
 */
export async function fetchSpliceAIPrediction(chrom, pos, ref, alt, onProgress = null) {
    const cleanChr = String(chrom).replace(/^chr/i, '').trim();
    const cleanRef = (ref || 'N').toUpperCase().trim();
    const cleanAlt = (alt || 'N').toUpperCase().trim();
    const endPos = pos + cleanRef.length - 1;

    if (onProgress) onProgress("Consultando predicción SpliceAI...");

    let ds_ag = 0, ds_al = 0, ds_dg = 0, ds_dl = 0;
    let dp_ag = 0, dp_al = 0, dp_dg = 0, dp_dl = 0;
    let found = false;
    let endpointUsed = "";

    // 1. Primary endpoint: Ensembl VEP REST API with SpliceAI plugin
    try {
        const vepUrl = `https://rest.ensembl.org/vep/human/region/${cleanChr}:${pos}:${endPos}:1/${cleanAlt}?SpliceAI=1`;
        const res = await fetch(vepUrl, {
            headers: { "Content-Type": "application/json", "Accept": "application/json" }
        });
        if (res.ok) {
            const data = await res.json();
            if (Array.isArray(data) && data.length > 0) {
                const tcList = data[0].transcript_consequences || [];
                for (const tc of tcList) {
                    if (tc.spliceai && typeof tc.spliceai === 'object') {
                        const sp = tc.spliceai;
                        ds_ag = parseFloat(sp.DS_AG ?? sp.ds_ag ?? 0);
                        ds_al = parseFloat(sp.DS_AL ?? sp.ds_al ?? 0);
                        ds_dg = parseFloat(sp.DS_DG ?? sp.ds_dg ?? 0);
                        ds_dl = parseFloat(sp.DS_DL ?? sp.ds_dl ?? 0);

                        dp_ag = parseInt(sp.DP_AG ?? sp.dp_ag ?? 0, 10);
                        dp_al = parseInt(sp.DP_AL ?? sp.dp_al ?? 0, 10);
                        dp_dg = parseInt(sp.DP_DG ?? sp.dp_dg ?? 0, 10);
                        dp_dl = parseInt(sp.DP_DL ?? sp.dp_dl ?? 0, 10);
                        found = true;
                        endpointUsed = "ensembl_vep";
                        break;
                    }
                }
            }
        }
    } catch (e) {
        console.warn("Ensembl VEP SpliceAI fetch warning:", e);
    }

    // 2. Fallback: GeneBE API
    if (!found) {
        let rawData = null;
        try {
            const shieldUrl = `${GENEBE_API_URL}/shield/v1/variant/hg38/${cleanChr}-${pos}-${cleanRef}-${cleanAlt}`;
            const res = await fetch(shieldUrl, {
                headers: { "Accept": "application/json" }
            });
            if (res.ok) {
                rawData = await res.json();
                endpointUsed = "shield";
            }
        } catch (e) {
            console.warn("GeneBE shield endpoint warning:", e);
        }

        if (!rawData) {
            try {
                const publicUrl = `${GENEBE_API_URL}/cloud/api-public/v1/variant?chr=${cleanChr}&pos=${pos}&ref=${cleanRef}&alt=${cleanAlt}&genome=hg38`;
                const res = await fetch(publicUrl, {
                    headers: { "Accept": "application/json" }
                });
                if (res.ok) {
                    const json = await res.json();
                    if (json && json.variants && json.variants.length > 0) {
                        rawData = json.variants[0];
                        endpointUsed = "public";
                    }
                }
            } catch (e) {
                console.warn("GeneBE public variant endpoint warning:", e);
            }
        }

        if (rawData) {
            if (rawData.spliceai && typeof rawData.spliceai === 'object') {
                ds_ag = parseFloat(rawData.spliceai.ds_ag ?? rawData.spliceai.DS_AG ?? 0);
                ds_al = parseFloat(rawData.spliceai.ds_al ?? rawData.spliceai.DS_AL ?? 0);
                ds_dg = parseFloat(rawData.spliceai.ds_dg ?? rawData.spliceai.DS_DG ?? 0);
                ds_dl = parseFloat(rawData.spliceai.ds_dl ?? rawData.spliceai.DS_DL ?? 0);

                dp_ag = parseInt(rawData.spliceai.dp_ag ?? rawData.spliceai.DP_AG ?? 0, 10);
                dp_al = parseInt(rawData.spliceai.dp_al ?? rawData.spliceai.DP_AL ?? 0, 10);
                dp_dg = parseInt(rawData.spliceai.dp_dg ?? rawData.spliceai.DP_DG ?? 0, 10);
                dp_dl = parseInt(rawData.spliceai.dp_dl ?? rawData.spliceai.DP_DL ?? 0, 10);
                found = true;
            } else if (rawData.ds_ag !== undefined || rawData.DS_AG !== undefined) {
                ds_ag = parseFloat(rawData.ds_ag ?? rawData.DS_AG ?? 0);
                ds_al = parseFloat(rawData.ds_al ?? rawData.DS_AL ?? 0);
                ds_dg = parseFloat(rawData.ds_dg ?? rawData.DS_DG ?? 0);
                ds_dl = parseFloat(rawData.ds_dl ?? rawData.DS_DL ?? 0);

                dp_ag = parseInt(rawData.dp_ag ?? rawData.DP_AG ?? 0, 10);
                dp_al = parseInt(rawData.dp_al ?? rawData.DP_AL ?? 0, 10);
                dp_dg = parseInt(rawData.dp_dg ?? rawData.DP_DG ?? 0, 10);
                dp_dl = parseInt(rawData.dp_dl ?? rawData.DP_DL ?? 0, 10);
                found = true;
            } else if (rawData.spliceai_max_score !== undefined || rawData.splice_score_selected !== undefined) {
                const maxVal = parseFloat(rawData.spliceai_max_score ?? rawData.splice_score_selected ?? 0);
                const effectStr = String(rawData.effect || rawData.splice_prediction_selected || '').toLowerCase();
                
                if (effectStr.includes('donor')) {
                    ds_dl = maxVal;
                    dp_dl = -1; // Donor junction is 1 bp upstream from +1
                } else if (effectStr.includes('acceptor')) {
                    ds_al = maxVal;
                    dp_al = -2; // Acceptor junction is 2 bp downstream from -2
                } else {
                    ds_dl = maxVal;
                    dp_dl = -1;
                }
                found = true;
            }
        }
    }

    if (!found) {
        return {
            success: false,
            available: false,
            message: "No se encontraron predicciones de SpliceAI para esta coordenada."
        };
    }

    // Safety checks for NaN
    ds_ag = Number.isFinite(ds_ag) ? Math.max(0, Math.min(1, ds_ag)) : 0;
    ds_al = Number.isFinite(ds_al) ? Math.max(0, Math.min(1, ds_al)) : 0;
    ds_dg = Number.isFinite(ds_dg) ? Math.max(0, Math.min(1, ds_dg)) : 0;
    ds_dl = Number.isFinite(ds_dl) ? Math.max(0, Math.min(1, ds_dl)) : 0;

    const maxScore = Math.max(ds_ag, ds_al, ds_dg, ds_dl);
    const maxLoss = Math.max(ds_al, ds_dl);
    const maxGain = Math.max(ds_ag, ds_dg);

    // Build recommendation connecting to Visual Cases
    let recommendation = {
        action: 'neutral',
        badge: 'Bajo Impacto',
        badgeClass: 'badge-neutral',
        title: 'Baja probabilidad de alteración de splicing canónico',
        text: 'Todos los puntajes delta de SpliceAI son inferiores al umbral de relevancia clínica (< 0.2). Es poco probable que la variante induzca cambios significativos en el patrón de empalme.',
        suggestedCase: null
    };

    if (maxScore >= 0.2) {
        if (maxLoss >= maxGain) {
            // Loss dominates -> Exon Skipping (Case 1) or Intron Retention (Case 3)
            const isVeryHigh = maxLoss >= 0.8;
            const isHigh = maxLoss >= 0.5;
            const impactLevel = isVeryHigh ? 'Muy Fuerte (> 0.8)' : (isHigh ? 'Fuerte (> 0.5)' : 'Sospechoso (> 0.2)');
            const lossType = ds_dl >= ds_al ? 'Donante (Donor Loss)' : 'Aceptor (Acceptor Loss)';

            recommendation = {
                action: 'loss',
                badge: `Pérdida ${lossType} (${maxLoss.toFixed(2)})`,
                badgeClass: isVeryHigh ? 'badge-danger' : (isHigh ? 'badge-warning' : 'badge-warning'),
                title: `Pérdida de Sitio de Splicing Canónico — Impacto ${impactLevel}`,
                text: `SpliceAI predice una pérdida ${lossType.toLowerCase()} significativa (DS: ${maxLoss.toFixed(2)}). El evento molecular más probable es el <strong>Salto de Exón (Exon Skipping - Caso 1)</strong> o alternativamente la <strong>Retención de Intrón (Caso 3)</strong>.`,
                suggestedCase: 1
            };
        } else {
            // Gain dominates -> Cryptic Splice Site (Case 2)
            const isVeryHigh = maxGain >= 0.8;
            const isHigh = maxGain >= 0.5;
            const impactLevel = isVeryHigh ? 'Muy Fuerte (> 0.8)' : (isHigh ? 'Fuerte (> 0.5)' : 'Sospechoso (> 0.2)');
            const gainType = ds_dg >= ds_ag ? 'Donante Críptico (Donor Gain)' : 'Aceptor Críptico (Acceptor Gain)';
            const deltaPos = ds_dg >= ds_ag ? dp_dg : dp_ag;
            const posSign = deltaPos > 0 ? `+${deltaPos}` : `${deltaPos}`;

            recommendation = {
                action: 'gain',
                badge: `Ganancia ${gainType} (${maxGain.toFixed(2)})`,
                badgeClass: isVeryHigh ? 'badge-danger' : (isHigh ? 'badge-warning' : 'badge-warning'),
                title: `Activación de Sitio Críptico — Impacto ${impactLevel}`,
                text: `SpliceAI predice la ganancia de un ${gainType.toLowerCase()} (DS: ${maxGain.toFixed(2)}). El evento molecular más probable es la <strong>Activación de Sitio Críptico (Caso 2)</strong>. La nueva unión ocurriría a <strong>${posSign} pb</strong> de la posición de la variante.`,
                suggestedCase: 2,
                deltaPos
            };
        }
    }

    return {
        success: true,
        available: true,
        endpointUsed,
        ds_ag, dp_ag,
        ds_al, dp_al,
        ds_dg, dp_dg,
        ds_dl, dp_dl,
        maxScore,
        recommendation
    };
}
