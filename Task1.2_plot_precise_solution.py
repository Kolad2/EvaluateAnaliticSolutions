import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow, c_tr_const
from scipy.optimize import curve_fit
from scipy.stats import lognorm
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def setup_custom_font(font_size=18):
	font_path = Path(".") / "assets" / "timesnewromanpsmt.ttf"

	if not font_path.exists():
		raise FileNotFoundError(f"Шрифт не найден: {font_path}")

	# Регистрируем шрифт в matplotlib
	mpl.font_manager.fontManager.addfont(str(font_path))

	# Получаем внутреннее имя семейства шрифта
	font_prop = mpl.font_manager.FontProperties(fname=str(font_path))
	font_name = font_prop.get_name()

	# Применяем глобально
	plt.rcParams.update({
		"font.family": font_name,
		"font.size": font_size,
		"axes.labelsize": font_size,
		"xtick.labelsize": font_size,
		"ytick.labelsize": font_size,
		"legend.fontsize": font_size,
		"figure.titlesize": font_size,

		# Для формул в $...$ — стиль близкий к Times
		"mathtext.fontset": "stix",
		"mathtext.rm": "STIXGeneral",
		"mathtext.it": "STIXGeneral:italic",
		"mathtext.bf": "STIXGeneral:bold",
	})

	return font_prop


def main():
	FONT_SIZE = 18
	setup_custom_font(FONT_SIZE)

	q = 0.90
	n_max = 60
	t = np.array([5.0, 25.0*q**(-16), 25.0*q**(-40)])
	lt = np.array([20, 50, 100])
	n = np.arange(0, n_max)

	c_q = [c_tr_pow(n_max, _t, q) for _t in t]
	c = [c_tr_const(n_max, 0.5*_t) for _t in lt]
	marker_type = ["o", "s", "D"]

	fig = plt.figure(figsize=(10, 5))
	ax = [fig.add_subplot(1, 2, i + 1) for i in range(2)]

	for i, _c in enumerate(c):
		_c = mask_keep_gap_edges(_c, 0.001)
		_fit_n, _fit_c, rmse = get_lognorm_fit(n, _c)
		ax[0].plot(-_fit_n, _fit_c, "-", linewidth=2, color="blue")
		ax[0].plot(-n, _c, marker_type[i], color="black", markersize=3)

	for i, _c in enumerate(c_q):
		_c = mask_keep_gap_edges(_c)
		print(np.nansum(_c))
		_fit_n, _fit_c, rmse = get_lognorm_fit(n, _c)
		print("error: ", rmse)
		ax[1].plot(-_fit_n, _fit_c, "-", linewidth=2, color="blue")
		ax[1].plot(-n, _c, marker_type[i], color="black", markersize=3)

	n_max_10 = np.ceil(n_max * np.log(2) / np.log(10))
	xticks_10 = np.arange(0, -n_max_10 - 1, -4)
	xticks = xticks_10 * (np.log(10) / np.log(2))

	ax[1].set_xlim(ax[0].get_xlim())
	ax[0].set_xticks(xticks)
	ax[1].set_xticks(xticks)

	xticklabels = [f'$10^{{{int(l)}}}$' for l in xticks_10]
	ax[0].set_xticklabels(xticklabels)
	ax[1].set_xticklabels(xticklabels)

	ax[0].set_xlabel("Relative size, $x/x_0$")
	ax[1].set_xlabel("Relative size, $x/x_0$")
	ax[0].set_ylabel("PDF")
	ax[1].set_ylabel("PDF")

	plt.tight_layout()
	plt.savefig("pictures/1.png", dpi=300)
	plt.show()


def mask_keep_gap_edges(y, thr=0.01):
	good = y > thr
	bad = ~good

	prev_good = np.r_[False, good[:-1]]
	next_good = np.r_[good[1:], False]

	bad_touching_good = bad & (prev_good | next_good)

	keep = good | bad_touching_good
	return np.where(keep, y, np.nan)


def lognorm_model(x, s, scale):
	return lognorm.pdf(x, s=s, loc=0, scale=scale)


def get_lognorm_fit(n, c, mask_func=None, thr=0.01, n_fit_points=300):
	n = np.asarray(n, dtype=float)
	c = np.asarray(c, dtype=float).copy()

	if mask_func is not None:
		c = mask_func(c, thr=thr)

	mask = np.isfinite(n) & np.isfinite(c) & (n > 0) & (c > 0)

	x_data = n[mask]
	y_data = c[mask]

	if len(x_data) < 3:
		return np.array([]), np.array([]), None

	s0 = 1.0
	scale0 = x_data[np.argmax(y_data)]

	try:
		popt, _ = curve_fit(
			lognorm_model,
			x_data,
			y_data,
			p0=[s0, scale0],
			bounds=([1e-8, 1e-8], [np.inf, np.inf]),
			maxfev=10000,
		)
	except RuntimeError:
		return np.array([]), np.array([]), None

	fit_n = np.linspace(x_data.min(), x_data.max(), n_fit_points)
	fit_c = lognorm_model(fit_n, *popt)
	residuals = y_data - lognorm_model(x_data, *popt)
	rmse = np.sqrt(np.mean(residuals ** 2))
	return fit_n, fit_c, rmse


if __name__ == "__main__":
	main()