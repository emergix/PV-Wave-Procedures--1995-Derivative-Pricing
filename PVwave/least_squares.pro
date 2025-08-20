common swaption_fit_common, periods, dates, tf, notionals, fixed_rates, type, $
       expiration, american, sliding, prices, zero, steps


pro test_swaptions
common swaption_fit_common, periods, dates, tf, notionals, fixed_rates, type, $
       expiration, american, sliding, prices, zero, steps
  periods = [6, 6, 6, 6]
  dates = [[0.5, 1., 1.5, 2., 2.5, 3.], $
           [0.5, 1., 1.5, 2., 2.5, 3.], $
           [0.5, 1., 1.5, 2., 2.5, 3.], $
           [0.5, 1., 1.5, 2., 2.5, 3.]]
  tf = [0.5, 0.5, 0.5, 0.5]
  notionals = [[1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0], $
               [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0], $
               [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0], $
               [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]]
  fixed_rates = [[0.12, 0.12, 0.12, 0.12, 0.12, 0.12], $
		 [0.12, 0.12, 0.12, 0.12, 0.12, 0.12], $
		 [0.12, 0.12, 0.12, 0.12, 0.12, 0.12], $
		 [0.12, 0.12, 0.12, 0.12, 0.12, 0.12]]
  type = [0, 0, 0, 0]
  expiration = [0.5, 1.0, 1.5, 2.0]
  american = [0, 0, 0, 0]
  sliding = [0, 0, 0, 0]
  prices = [14.63, 11.22, 7.74, 4.42]
  zero = transpose ([[1., 0.1], [2., 0.105], [3., 0.11], [4., 0.1125], $
                     [5, 0.115], [6., 0.117], [7., 0.123]])
end

pro init_swaptions, zero_curve, vb_swaptions
common swaption_fit_common, periods, dates, tf, notionals, fixed_rates, type, $
       expiration, american, sliding, prices, zero, steps
  sa = size (vb_swaptions)
  s = sa (1)
  periods = fltarr (s)
  dates = fltarr (12, s)
  tf = fltarr (s)
  notionals = fltarr (12, s)
  fixed_rates = fltarr (12, s)
  type = fltarr (s)
  expiration = fltarr (s)
  american = fltarr (s)
  sliding = fltarr (s)
  prices = fltarr (s)
  
  for i = 0, s - 1 do begin
    n = fix (vb_swaptions (i, 3) / vb_swaptions (i, 2))
    periods (i) = n
    for j = 0, n - 1 do begin
      dates (j, i) = (j + 1) * (vb_swaptions (i, 2) / 12)
      notionals (j, i) = 1.0
      fixed_rates (j, i) = vb_swaptions (i, 1)
      tf (i) = 0.5
      type (i) = vb_swaptions (i, 0)
      expiration (i) = vb_swaptions (i, 4)
      american (i) = vb_swaptions (i, 5)
      sliding (i) = vb_swaptions (i, 6)
      prices (i) = vb_swaptions (i, 7)
    end
  end
  zero = zero_curve
end

function swaption_fit, m, x
common swaption_fit_common, periods, dates, tf, notionals, fixed_rates, type, $
       expiration, american, sliding, prices, zero, steps
  n = 15
  sigma = x (0)
  a = x (1)
  result = fltarr (m)
  if (sigma GT 0.15) or (sigma LT 0.005) or (abs (a) GT 0.1) then begin
    result (*) = 1000
    return,result
  end
  vb_print, 'sigma = ' + string (sigma) + ', a = ' + string (a)
  crbarbre = HWbuild_crbtree (zero, a, sigma, n)
  for i = 1, m - 1 do begin
    hw_price = HWSwapOption (crbarbre, periods (i), dates (*, i), tf (i), $
                             notionals (*, i), fixed_rates (*, i), type (i), $
                             expiration (i), american (i), $
                             sliding (i))
    result (i) = hw_price - prices (i)
  end
  result = 10 * result
  print, result
  return, result
end

function fit_interest_parameters, n
common swaption_fit_common, periods, dates, tf, notionals, fixed_rates, type, $
       expiration, american, sliding, prices, zero, steps
  steps = n
  r = nlinlsq ("swaption_fit", n_elements (periods), 2, $
               XGuess = [0.02, 0.01])
  sigma = r (0)
  a = r (1)
  n = steps
  crbarbre = HWbuild_crbtree (zero, a, sigma, n)
  return, {, a:a, sigma:sigma, theta:crbarbre.arbre.theta}
end

pro test_price
common swaption_fit_common, periods, dates, tf, notionals, fixed_rates, type, $
       expiration, american, sliding, prices, zero, steps
 a = 0.03
 sigma = 0.015
;  a = 0.0099
;  sigma = 0.019
  n = steps
  crbarbre = HWbuild_crbtree (zero, a, sigma, n)
  for i = 0, n_elements (periods) - 1 do begin
    print, 'prix swaption= ', $
      HWSwapOption (crbarbre, periods (i), dates (*, i), tf (i), $
                    notionals (*, i), fixed_rates (*, i), $
	            type (i), $
                    expiration (i), american (i), $
                    sliding (i))
  end
end

function index_options_fit, ops, ind, zerocurve, yields, dates, steps, aa, ssigma
common option_fit_common, yieldlist, t_divlist, a, sigma, options, n_options, zero, index, crbtree
  index = ind
  yieldlist = yields
  t_divlist = dates
  a = aa
  sigma = ssigma
  ss = size (ops)
  s = ss (1)
  options = fltarr (80, s)
  for i = 0, s - 1 do begin
    options (1, i) = 600
    if ops (i, 1) eq 1 then options (1, i) = options (1, i) + 50
    if ops (i, 5) eq 1 then options (1, i) = options (1, i) + 1
    if ops (i, 2) eq 1 then options (1, i) = options (1, i) + 2    
    options (2, i) = ops (i, 1)
    options (3, i) = ops (i, 3)
    if ops (i, 1) eq 1 then options (4, i) = 0 else options (4, i) = ops (i, 3)
    options (10, i) = ops (i, 1)
    options (11, i) = ops (i, 2)
  end
  n_options = s
  zero = zerocurve
  nb = 8
  crbtree = HWbuild_crbtree (zero, a, sigma, nb)

  for i = 0, n_options - 1 do begin
    print, 'option ', i
    pm, options (*, i)
  end

  vol = 0.2
  corr = +0.0

  hwlnprob = HWLNtree (crbtree, vol, sigma, corr, index)
  result = HWLNOption (crbtree, hwlnprob, options, n_options, $
                       yieldlist, t_divlist)

  print, result

end