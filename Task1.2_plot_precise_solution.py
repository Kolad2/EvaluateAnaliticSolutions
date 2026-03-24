import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from qmath import QBinom, QPochhammer
from solutions import c_tr_pow, c_tr_const
from scipy.optimize import curve_fit
from scipy.stats import lognorm


def main():
	q = 0.90
	t = 0
	n_max = 60
	t = np.array([1.0, 5.0, 25.0*q**(-30), 25.0*q**(-45)])
	lt = np.array([20, 50, 100])
	n = np.arange(0, n_max)
	#%%
	c_q = [c_tr_pow(n_max, _t, q) for _t in t]
	c = [c_tr_const(n_max, 0.5*_t) for _t in lt]

	fig = plt.figure(figsize=(10, 5))
	ax = [fig.add_subplot(1,2,i+1) for i in range(2)]
	for _c in c:
		_c = mask_keep_gap_edges(_c, 0.001)
		_fit_n, _fit_c, rmse = get_lognorm_fit(n, _c)
		#print(np.nansum(_c))
		ax[0].plot(_fit_n, _fit_c, "-", linewidth=2, color="blue")
		ax[0].plot(n, _c, "-o", color="black", markersize=3)
	for _c in c_q:
		_c = mask_keep_gap_edges(_c)
		print(np.nansum(_c))
		_fit_n, _fit_c, rmse = get_lognorm_fit(n, _c)
		print("error: ", rmse)
		ax[1].plot(_fit_n, _fit_c, "-", linewidth=2, color="blue")
		ax[1].plot(n, _c, "-o", color="black", markersize=3)

	n_max_10 = np.ceil(n_max * np.log(2) / np.log(10))
# xticks_10 = np.arange(0, -n_max_10-1, -4)
	# xticks = xticks_10*(np.log(10)/np.log(2))
	# ax[0].set_xticks(xticks)
	# ax[1].set_xticks(xticks)
	# xticklabels = [f'$10^{{{int(l)}}}$' for l in xticks_10]
	# ax[0].set_xticklabels(xticklabels)
	# ax[1].set_xticklabels(xticklabels)
	#
	n_max_10 = np.ceil(n_max * np.log(2) / np.log(10))
	# yticks_10 = np.arange(0, n_max_10 - 1, 4)
	# yticks = yticks_10 * (np.log(10) / np.log(2))
	# ax[0].set_yticks(yticks)
	# ax[1].set_yticks(yticks)
	# yticklabels = [f'$10^{{{int(l)}}}$' for l in yticks_10]
	# ax[0].set_yticklabels(yticklabels)
	# ax[1].set_yticklabels(yticklabels)
	#
	ax[0].set_xlabel("Относительный размер, $x/x_0$")
	ax[1].set_xlabel("Относительный размер, $x/x_0$")
	ax[0].set_ylabel("Концентрация частиц, $c(x,t)/N_0$")
	ax[1].set_ylabel("Концентрация частиц, $c(x,t)/N_0$")
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
		popt, _ = curve_fit(lognorm_model, x_data, y_data, p0=[s0, scale0], bounds=([1e-8, 1e-8], [np.inf, np.inf]), maxfev=10000,)
	except RuntimeError:
		return np.array([]), np.array([]), None

	fit_n = np.linspace(x_data.min(), x_data.max(), n_fit_points)
	fit_c = lognorm_model(fit_n, *popt)
	residuals = y_data - lognorm_model(x_data, *popt)
	rmse = np.sqrt(np.mean(residuals ** 2))
	return fit_n, fit_c, rmse

if __name__ == "__main__":
	main()