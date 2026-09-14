"""Agrega resultados auditados y dibuja datos frente a predicciones sin ajuste.

No calcula una significancia de descubrimiento: aún faltan sistemáticas y una
estimación validada de fondos instrumentales. El chi-cuadrado es un diagnóstico
condicional del acuerdo estadístico, no una prueba de entrelazamiento.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2

from analysis import BASE, RESULTS, BINS, add_hist, normalized_shape, pipeline_hashes

LABELS = {'data': 'Colisiones reales', 'ttbar': r'$t\bar{t}$', 'single_top': 'Top individual',
          'diboson': 'Dos bosones', 'Zjets': 'Z + jets', 'Wjets': 'W + jets',
          'ttV': r'$t\bar{t}$ + bosones', 'HWW': r'Higgs $\to WW$'}
COLORS = {'ttbar': '#4b83b5', 'single_top': '#e6a657', 'diboson': '#70a78a',
          'Zjets': '#bd83b8', 'Wjets': '#d17474', 'ttV': '#bdc766', 'HWW': '#929292'}
ORDER = ['HWW', 'ttV', 'Wjets', 'Zjets', 'diboson', 'single_top', 'ttbar']
RUN_LABEL = ''
XLABELS = {'electron_pt': r'$p_T(e)$ [GeV]', 'muon_pt': r'$p_T(\mu)$ [GeV]',
           'n_jets': 'Número de jets aceptados', 'delta_phi': r'$|\Delta\phi(e,\mu)|$ [rad]',
           'delta_eta': r'$|\Delta\eta(e,\mu)|$ (sin unidades)',
           'm_emu': r'$m(e,\mu)$ [GeV]', 'met': 'Momento transversal faltante [GeV]',
           'n_bjets': 'Número de jets con etiqueta b'}


def poisson_errors(counts):
    """Intervalos de Garwood al 68.27% para conteos Poisson; válidos en cero."""
    n = np.asarray(counts, float)
    alpha = 1 - .682689492137
    lower = np.zeros_like(n)
    positive = n > 0
    lower[positive] = .5 * chi2.ppf(alpha / 2, 2 * n[positive])
    upper = .5 * chi2.ppf(1 - alpha / 2, 2 * (n + 1))
    return np.array([n - lower, upper - n])


def aggregate(run):
    """Suma procesos disjuntos y conserva los archivos individuales en JSON."""
    groups = {}
    for path in run['files']:
        sample = json.loads((BASE / path).read_text())
        assert sample['code_sha256'] == run['code_sha256'], 'Se intentó mezclar versiones de código'
        assert sample['config'] == run['config'], 'Se intentó mezclar selecciones'
        group = groups.setdefault(sample['group'], {'entries': 0, 'files': 0, 'cutflow': {},
                                                    'regions': {}, 'histograms': {}, 'negative_weights': 0})
        group['entries'] += sample['entries']
        group['files'] += 1
        group['negative_weights'] += sample['negative_final_weights']
        for category in ('cutflow', 'regions'):
            for label, stat in sample[category].items():
                dest = group[category].setdefault(label, {'raw': 0, 'sumw': 0., 'sumw2': 0.})
                for key in dest:
                    dest[key] += stat[key]
        for region, hists in sample['histograms'].items():
            dest = group['histograms'].setdefault(region, {})
            for obs, hist in hists.items():
                if obs in dest:
                    add_hist(dest[obs], hist)
                else:
                    dest[obs] = hist
    return groups


def total_mc(groups, region, variable):
    """Devuelve la suma nominal y su varianza MC, conservando pesos firmados."""
    n = len(BINS[variable]) - 1
    h, v = np.zeros(n), np.zeros(n)
    for g in ORDER:
        if g in groups:
            hist = groups[g]['histograms'][region][variable]
            h += hist['sumw']
            v += hist['sumw2']
    return h, v


def comparison_panel(ax, ratio, groups, region, variable):
    """Gráfico de conteos absolutos; la banda MC y los errores de datos son distintos."""
    bins = BINS[variable]
    centres = (bins[:-1] + bins[1:]) / 2
    bottom = np.zeros(len(centres))
    for g in ORDER:
        if g not in groups:
            continue
        h = np.array(groups[g]['histograms'][region][variable]['sumw'])
        # No se recortan bins negativos para hacer más bonito el gráfico.
        ax.stairs(bottom+h, bins, baseline=bottom, fill=True, color=COLORS[g], label=LABELS[g], linewidth=.4)
        bottom += h
    mc, var = total_mc(groups, region, variable)
    err = np.sqrt(var)
    ax.stairs(mc+err, bins, baseline=mc-err, fill=True, facecolor='none',
              edgecolor='#444444', hatch='////', linewidth=0, label='Estadística MC')
    data = np.array(groups['data']['histograms'][region][variable]['raw'])
    de = poisson_errors(data)
    ax.errorbar(centres, data, yerr=de, fmt='o', color='#161616', markersize=3, linewidth=.8, label='Datos reales')
    ax.set_ylabel('Eventos / intervalo')
    ax.set_ylim(bottom=0, top=max(float(np.max(data+de[1])), float(np.max(mc+err)), 1)*1.28)
    ax.grid(axis='y', alpha=.14)
    valid = mc > 0
    relative_error = np.divide(err, mc, out=np.zeros_like(err), where=valid)
    ratio.stairs(1+relative_error, bins, baseline=1-relative_error, fill=True, facecolor='#cccccc', linewidth=0)
    ratio.errorbar(centres[valid], data[valid]/mc[valid], yerr=de[:, valid]/mc[valid], fmt='o', color='black', markersize=3, linewidth=.8)
    ratio.axhline(1, color='#4b5563', lw=.8)
    ratio.set_ylabel('Datos / MC')
    ratio.set_xlabel(XLABELS[variable])
    if valid.any():
        lower = min(.8, float(np.min((data[valid]-de[0, valid])/mc[valid])) - .04,
                    float(np.min(1-relative_error[valid])) - .04)
        upper = max(1.2, float(np.max((data[valid]+de[1, valid])/mc[valid])) + .04,
                    float(np.max(1+relative_error[valid])) + .04)
        ratio.set_ylim(max(0, lower), upper)
    else:
        ratio.set_ylim(0, 2)
    ratio.grid(axis='y', alpha=.2)
    if variable in ('n_jets', 'n_bjets'):
        ticks = centres.astype(int)
        ratio.set_xticks(ticks, [str(x) for x in ticks[:-1]] + [f'>={ticks[-1]}'])
    elif variable not in ('delta_phi', 'delta_eta'):
        ratio.text(.98, -.62, f'Último intervalo incluye valores >= {bins[-1]:g}',
                   transform=ratio.transAxes, ha='right', fontsize=7)


def draw_comparisons(groups, region, out):
    fig = plt.figure(figsize=(12, 9.5), layout='constrained')
    grid = fig.add_gridspec(4, 2, height_ratios=[3, 1, 3, 1])
    axes = []
    for i, variable in enumerate(('electron_pt', 'muon_pt', 'n_jets', 'delta_phi')):
        row, col = 2*(i//2), i % 2
        ax = fig.add_subplot(grid[row, col])
        ratio = fig.add_subplot(grid[row+1, col], sharex=ax)
        plt.setp(ax.get_xticklabels(), visible=False)
        comparison_panel(ax, ratio, groups, region, variable)
        axes.append(ax)
    axes[0].legend(fontsize=9, ncols=3, loc='upper right')
    tag = 'al menos una' if region == 'SR1b' else 'al menos dos'
    fig.suptitle(RUN_LABEL + f'ATLAS Open Data: e-mu, cargas opuestas, {tag} etiqueta(s) b\n'
                 '2015-2016 | 13 TeV | L = 36 fb$^{-1}$ (redondeada) | predicción nominal sin ajuste', fontsize=12)
    fig.savefig(out / f'datos_mc_{region}.png', dpi=175)
    plt.close(fig)


def draw_angles(groups, region, out, normalized=False):
    """Dos separaciones leptónicas, con incertidumbres sólo estadísticas.

    Conteos: panel inferior datos/MC. Forma normalizada: fracción por intervalo;
    la incertidumbre diagonal procede de la covarianza completa de normalización.
    Normalizar elimina la información del rendimiento total, no los fondos.
    """
    fig = plt.figure(figsize=(11, 5.8), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[3, 1])
    for col, obs in enumerate(('delta_phi', 'delta_eta')):
        ax = fig.add_subplot(grid[0, col])
        lower = fig.add_subplot(grid[1, col], sharex=ax)
        if not normalized:
            comparison_panel(ax, lower, groups, region, obs)
        else:
            edges = BINS[obs]
            centres = (edges[:-1] + edges[1:]) / 2
            data = np.asarray(groups['data']['histograms'][region][obs]['raw'], float)
            mc, var = total_mc(groups, region, obs)
            pd, cd = normalized_shape(data, data)
            pm, cm = normalized_shape(mc, var)
            ed, em = np.sqrt(np.maximum(np.diag(cd), 0)), np.sqrt(np.maximum(np.diag(cm), 0))
            ax.stairs(pm, edges, color='#3476a8', label='MC nominal + fondos')
            ax.stairs(pm+em, edges, baseline=pm-em, fill=True, color='#3476a8', alpha=.2,
                      label='Incertidumbre estadística MC')
            ax.errorbar(centres, pd, yerr=ed, fmt='o', color='black', markersize=4,
                        label='Colisiones reales')
            ax.set_ylabel('Fracción de eventos / intervalo')
            valid = pm > 0
            lower.errorbar(centres[valid], pd[valid]/pm[valid], yerr=ed[valid]/pm[valid],
                           fmt='o', color='black', markersize=3)
            band = np.divide(em, pm, out=np.zeros_like(pm), where=valid)
            lower.stairs(1+band, edges, baseline=1-band, fill=True, color='#3476a8', alpha=.2)
            lower.axhline(1, color='gray', lw=.8)
            lower.set_ylabel('Datos / MC')
            lower.set_xlabel(XLABELS[obs])
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.grid(axis='y', alpha=.15)
        lower.grid(axis='y', alpha=.15)
        if col == 0:
            handles, labels = ax.get_legend_handles_labels()
    tag = 'al menos 1' if region == 'SR1b' else 'al menos 2'
    mode = 'formas normalizadas' if normalized else 'conteos absolutos'
    fig.suptitle(RUN_LABEL + f'ATLAS Open Data | e-mu, {tag} b-tags | {mode}\n'
                 '2015-2016 | 13 TeV | 36 fb$^{-1}$ (redondeada) | nivel reconstruido', fontsize=11)
    fig.legend(handles, labels, loc='outside lower center', ncols=3, fontsize=8)
    stem = f'angulos_{"formas" if normalized else "conteos"}_{region}'
    fig.savefig(out / f'{stem}.png', dpi=175)
    fig.savefig(out / f'{stem}.pdf')
    plt.close(fig)


def diagnostics(groups):
    """Evalúa rendimientos y formas. Todo el contraste es solo estadístico."""
    results = {}
    for region in ('PRETAG', 'ZEROb', 'SR1b', 'SR2b', 'SS1b'):
        ndata = groups['data']['regions'][region]['raw']
        mc = sum(g['regions'][region]['sumw'] for k, g in groups.items() if k != 'data')
        var = sum(g['regions'][region]['sumw2'] for k, g in groups.items() if k != 'data')
        signal = groups['ttbar']['regions'][region]['sumw']
        results[region] = {'data': ndata, 'mc': mc, 'mc_stat': float(np.sqrt(var)),
                           'ttbar_fraction_in_model': signal/mc if mc > 0 else None,
                           'data_over_mc': ndata/mc if mc > 0 else None, 'shapes': {}}
        for observable in BINS:
            data = np.array(groups['data']['histograms'][region][observable]['raw'], float)
            pred, v = total_mc(groups, region, observable)
            if pred.sum() <= 0 or data.sum() <= 0:
                continue
            p_data, c_data = normalized_shape(data, data)
            p_mc, c_mc = normalized_shape(pred, v)
            covariance = c_data + c_mc
            rank = int(np.linalg.matrix_rank(covariance))
            delta = p_data - p_mc
            value = float(delta @ np.linalg.pinv(covariance) @ delta)
            results[region]['shapes'][observable] = {'chi2_stat_shape': value, 'rank': rank,
                'p_approx_stat_only': float(chi2.sf(value, rank)) if rank else None,
                'covariance_data': c_data.tolist(), 'covariance_mc': c_mc.tolist(),
                'warning': 'Diagnóstico asintótico solo estadístico; sin sistemáticas, no es significancia de descubrimiento.'}
    return results


def draw_control(groups, out):
    fig = plt.figure(figsize=(11, 5), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[3, 1])
    for col, (region, obs) in enumerate((('PRETAG', 'n_bjets'), ('SS1b', 'delta_phi'))):
        ax, ratio = fig.add_subplot(grid[0, col]), fig.add_subplot(grid[1, col])
        comparison_panel(ax, ratio, groups, region, obs)
        ax.set_title('Antes de exigir etiqueta b' if region == 'PRETAG' else 'Control: cargas iguales', fontsize=11)
        if col == 0:
            handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8.5, ncols=5, loc='outside lower center')
    fig.suptitle(RUN_LABEL + 'Controles de selección y de modelado | incertidumbres estadísticas', fontsize=12)
    fig.savefig(out / 'controles.png', dpi=175)
    plt.close(fig)


def draw_extra(groups, out):
    """Controles adicionales que no se usan para retocar la selección."""
    fig = plt.figure(figsize=(11, 5), layout='constrained')
    grid = fig.add_gridspec(2, 2, height_ratios=[3, 1])
    for col, obs in enumerate(('m_emu', 'met')):
        ax, ratio = fig.add_subplot(grid[0, col]), fig.add_subplot(grid[1, col])
        comparison_panel(ax, ratio, groups, 'SR1b', obs)
        if col == 0:
            handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=8.5, ncols=5, loc='outside lower center')
    fig.suptitle(RUN_LABEL + 'SR1b: masa electrón-muón y momento faltante | errores estadísticos', fontsize=12)
    fig.savefig(out / 'variables_adicionales.png', dpi=175)
    plt.close(fig)


def main():
    global RUN_LABEL
    parser = argparse.ArgumentParser()
    parser.add_argument('--partial', action='store_true')
    args = parser.parse_args()
    RUN_LABEL = 'CORRIDA PARCIAL - faltan muestras\n' if args.partial else ''
    run = json.loads((RESULTS / ('run_partial.json' if args.partial else 'run.json')).read_text())
    assert run['code_sha256'] == hashlib.sha256((BASE / 'analysis.py').read_bytes()).hexdigest(), 'El código actual difiere del utilizado; vuelva a ejecutar el análisis'
    if not args.partial:
        assert run['complete'], 'No se pueden publicar figuras completas con una corrida parcial'
        assert run['global_checks']['all_checks_passed'], 'La corrida no completó los controles globales'
    assert run.get('pipeline_hashes') == pipeline_hashes(), 'Los controles importados cambiaron desde la corrida'
    groups = aggregate(run)
    summary = {'complete': run['complete'], 'config': run['config'], 'groups': groups,
               'diagnostics': diagnostics(groups), 'n_files': len(run['files'])}
    out = BASE / 'figures' / ('angular_v2_partial' if args.partial else 'angular_v2')
    out.mkdir(parents=True, exist_ok=True)
    for region in ('SR1b', 'SR2b'):
        draw_comparisons(groups, region, out)
        draw_angles(groups, region, out)
        draw_angles(groups, region, out, normalized=True)
    draw_control(groups, out)
    draw_extra(groups, out)
    (RESULTS / ('summary_partial.json' if args.partial else 'summary.json')).write_text(json.dumps(summary, indent=2), encoding='utf-8')
    with (RESULTS / 'cutflow.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['group', 'stage', 'raw', 'sumw', 'sumw2', 'stat_error'])
        for name, group in groups.items():
            for stage, stat in group['cutflow'].items():
                writer.writerow([name, stage, stat['raw'], stat['sumw'], stat['sumw2'], np.sqrt(stat['sumw2'])])
    # Tabla plana: permite comprobar cada barra sin depender de Matplotlib.
    with (RESULTS / ('histograms_partial.csv' if args.partial else 'histograms.csv')).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.writer(stream)
        writer.writerow(['group', 'region', 'observable', 'bin_low', 'bin_high',
                         'last_bin_includes_overflow', 'raw', 'sumw', 'sumw2'])
        for name, group in groups.items():
            for region, hists in group['histograms'].items():
                for variable, h in hists.items():
                    edges = BINS[variable]
                    for i in range(len(edges)-1):
                        overflow = i == len(edges)-2 and variable not in ('delta_phi', 'delta_eta')
                        writer.writerow([name, region, variable, edges[i], edges[i+1],
                                         overflow, h['raw'][i], h['sumw'][i], h['sumw2'][i]])
    print(json.dumps({r: {k: v for k,v in d.items() if k != 'shapes'} for r,d in summary['diagnostics'].items()}, indent=2))


if __name__ == '__main__':
    main()
