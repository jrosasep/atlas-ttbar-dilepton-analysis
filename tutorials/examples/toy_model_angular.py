"""Toy model de los observables angulares del canal electrón-muón.

Utiliza cinco eventos inventados para testear cómo se calcularan los observables angulares del sistema:
    1. La separación azimutal periódica (Delta phi).
    2. La separación absoluta en pseudorrapidez (|Delta eta|).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def observables_angulares(eta_e, phi_e, eta_mu, phi_mu):
    """Calcula Delta phi y |Delta eta| para uno o varios eventos.

    Parámetros
    ----------
    eta_e, eta_mu : números o arreglos de NumPy
        Pseudorrapideces reconstruidas del electrón y del muón.
    phi_e, phi_mu : números o arreglos de NumPy
        Ángulos azimutales del electrón y del muón, expresados en radianes.

    Retorna
    -------
    delta_phi : número o arreglo de NumPy
        Separación azimutal mínima. Sus valores quedan entre 0 y pi.
    delta_eta : número o arreglo de NumPy
        Valor absoluto de la diferencia de pseudorrapidez.
    """

    # la separación entre las particulas respectro a sus angulos polares, 
    # la información angular se compacta en la pseudorrapidez, finalmente|Delta eta|. 
    # No se necesita una corrección periódica puesto que thetha no es un
    # ángulo periódico. Sólo restamos ambos valores de la pseudorrapidez del electron y el muon
    # y tomamos el valor absoluto.
    delta_eta = np.abs(eta_e - eta_mu)

    # Esta diferencia puede producir valores fuera del intervalo [-pi, pi].
    # Por ejemplo, 179 grados menos -179 grados entrega 358 grados, aunque
    # ambas direcciones están físicamente separadas por sólo 2 grados.
    diferencia_phi = phi_e - phi_mu

    # atan2(sin(delta), cos(delta)) devuelve el ángulo equivalente dentro de
    # (-pi, pi]. El valor absoluto conserva únicamente la separación mínima.
    delta_phi = np.abs(
        np.arctan2(
            np.sin(diferencia_phi),
            np.cos(diferencia_phi),
        )
    )

    return delta_phi, delta_eta


def construir_eventos_toy():
    """Construye cinco eventos inventados con situaciones "sencillas"."""

    # Cada posición de los cuatro arreglos describe un mismo evento.
    # Por ejemplo, los cuatro valores de índice 0 pertenecen al evento 1.
    eta_e = np.array([0.2, 1.1, -0.5, 2.2, -1.0])
    eta_mu = np.array([0.3, -0.2, -0.5, 1.7, 0.8])

    # Es más sencillo inventar y leer los ángulos en grados.
    phi_e_grados = np.array([179.0, 10.0, 0.0, 90.0, -45.0])
    phi_mu_grados = np.array([-179.0, 20.0, 180.0, 90.0, 135.0])

    # NumPy evalúa seno, coseno y atan2 en radianes. Por eso, si usamos grados sexagesimales, 
    # convertimos los valores a radianes antes de calcular los observables, esto utilizando la funcion np.deg2rad.
    phi_e = np.deg2rad(phi_e_grados)
    phi_mu = np.deg2rad(phi_mu_grados)

    return eta_e, phi_e, eta_mu, phi_mu


def imprimir_eventos(eta_e, phi_e, eta_mu, phi_mu, delta_phi, delta_eta):
    """Muestra las entradas y los observables calculados para cada evento."""

    encabezado = (
        "Evento | eta_e | eta_mu | phi_e [grados] | phi_mu [grados] "
        "| Delta phi [grados] | |Delta eta|"
    )
    print(encabezado)
    print("-" * len(encabezado))

    for i in range(len(eta_e)):
        print(
            f"{i + 1:>6} | "
            f"{eta_e[i]:>5.2f} | "
            f"{eta_mu[i]:>6.2f} | "
            f"{np.rad2deg(phi_e[i]):>15.1f} | "
            f"{np.rad2deg(phi_mu[i]):>16.1f} | "
            f"{np.rad2deg(delta_phi[i]):>18.1f} | "
            f"{delta_eta[i]:>11.2f}"
        )


def graficar_histogramas(delta_phi, delta_eta, ruta_salida):
    """Construye y guarda los histogramas de los cinco eventos toy."""

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Histograma de Delta phi. Sus límites físicos son 0 y pi radianes.
    axes[0].hist(
        delta_phi,
        bins=np.linspace(0, np.pi, 7),
        color="royalblue",
        edgecolor="black",
    )
    axes[0].set_xlim(0, np.pi)
    axes[0].set_xlabel(r"$\Delta\phi_{e\mu}$ [rad]")
    axes[0].set_ylabel("Eventos")
    axes[0].set_title("Separación azimutal")

    # Histograma de |Delta eta|. Esta variable es adimensional.
    axes[1].hist(
        delta_eta,
        bins=np.linspace(0, 3, 7),
        color="darkorange",
        edgecolor="black",
    )
    axes[1].set_xlim(0, 3)
    axes[1].set_xlabel(r"$|\Delta\eta_{e\mu}|$")
    axes[1].set_ylabel("Eventos")
    axes[1].set_title("Separación en pseudorrapidez")

    fig.suptitle(r"toy model: observables angulares del sistema $e\mu$")
    fig.tight_layout()
    fig.savefig(ruta_salida, dpi=200, bbox_inches="tight")
    plt.show()


def main():
    """Ejecuta, en orden, todas las etapas del toy model."""


#137#
    # 1. Construir los eventos inventados.
    eta_e, phi_e, eta_mu, phi_mu = construir_eventos_toy()

    # Todos los arreglos deben contener el mismo número de eventos.
    if not (len(eta_e) == len(phi_e) == len(eta_mu) == len(phi_mu)):
        raise ValueError("Los arreglos de entrada deben tener la misma longitud.")

    # 2. Calcular un par de observables por evento.
    delta_phi, delta_eta = observables_angulares(
        eta_e,
        phi_e,
        eta_mu,
        phi_mu,
    )

    # 3. Verificar propiedades básicas que siempre deben cumplirse.
    assert np.all((0 <= delta_phi) & (delta_phi <= np.pi))
    assert np.all(delta_eta >= 0)

    # El evento 1 fue diseñado para verificar la periodicidad: 179 y -179
    # grados deben estar separados por 2 grados, no por 358 grados.
    assert np.isclose(np.rad2deg(delta_phi[0]), 2.0)

    # 4. Mostrar los valores de cada evento para poder revisarlos a mano.
    imprimir_eventos(
        eta_e,
        phi_e,
        eta_mu,
        phi_mu,
        delta_phi,
        delta_eta,
    )

    # 5. Guardar la figura junto al archivo Python.
    ruta_salida = Path(__file__).with_name("toy_model_histograms.png")
    graficar_histogramas(delta_phi, delta_eta, ruta_salida)
    print(f"\nFigura guardada en: {ruta_salida}")
main()
