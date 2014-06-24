PRO Landsat_invert_from_Aaron

; Weight Radiances
file_name = 'C:\Users\ejipci\Documents\Emmett_Files\Work\Landsat TIRS\Data From Aaron\Merged_C_400.txt'
;datafile  = READ_ASCII(file_name, DATA_START=0)
n_plus_matt = 6180    ; Ind. observations.  The first col. are Matt's wts.  need to skip them.
r           = 82      ; Prdictor variables
n           = n_plus_matt - 1

; Read in the Weight Radiances
data = FLTARR(n_plus_matt, r)
openr, 1, file_name
readf, 1, data
close, 1

; Plot some stuff
; Matts_wts = data[0,*]
;window, 0
;plot, data(0,*), title='Matts Weights'
;window, 1
;plot, data(5000,*), title='A Sample Set of Wt. Radiances'
setcolor
; Omit the first column of Ghost "weights"
L = FLTARR(n, r)
L = data[1:n, *]





; make L, n x (r+1).....r prdictor variables.  Pg. 378-79.
; n = 6179
; r = 82
; this is 6179 x 83 now...then add a first row of 'ones'....70 ones.
ones = fltarr(n)
ones = ones + 1.0  ;  make a vector of ones
L_ones = fltarr(n,(r+1))  ; put the vecdtor into L
L_ones[*,0] = ones

L_ones[*,1:r] = L   ; put the rest of L into the new L matrix


; ******************************************
; Bias Radiances
; Read in the L8_bias goes Res file  
; Col 0 = Latitude
; Col 1 = Longitude
; Col 2 = BIAS VALUES
; Col 3 = Pixel position
file_name = 'C:\Users\ejipci\Documents\Emmett_Files\Work\Landsat TIRS\Data From Aaron\Biases.txt'
datafile  = READ_ASCII(file_name, DATA_START=0)  ; 3 x 90
Bias = datafile.(0) *(-1)
;plot, Bias, title='Bias Radiances'


; Just make a smaller matrix so the SVD doesn't puke
k = 900
L_ones = L_ones[0:k-1,*]
Bias = Bias[0:k-1]
 



; Compute least square inverse using IDL  INVERT FUNCTION
;Lpsu1  = INVERT((transpose(L) ## L), /DOUBLE) ## transpose(L) 
;w_est1 = B ## Lpsu1
;B_est1 = transpose(L) ## w_est1


; Compute least squares inverse using SVD for inverse part w/ conditioning
print, '...Make square'
Lsquare = transpose(L_ones) ## L_ones  ; n x n
print, 'done.'
print, '...Decompose uisng SVD'
SVDC, Lsquare, w, U, V, /double
print, 'done.'


; make diag matrix.  Only elements with absolute values greater than 1E-5 are reciprocated
N = N_elements(w)
WP = fltarr(N,N)
FOR k=0, N-1 DO $
    IF abs(W(k)) GE 1.0E-2 THEN WP(k,k) = 1.0/W(k)
L_inv = V ## WP ## transpose(U)

Lpsu2  = L_inv ## transpose(L_ones)
beta = Bias ## Lpsu2
Bias_est = transpose(L_ones) ## beta
window, 4
plot, beta, title='Regression Coeff'


; Compare the two methods
window, 5
plot,  Bias, thick=1, title='Bias (Blk) vs Bias_est (Grn)';, yrange=[-1.1,-0.3], title
oplot, Bias_est, color=3, thick=1
;oplot, B_est1, color=2, thick=1
;oplot, SMOOTH(B_est2,2), color=4, thick=2
;oplot, B, thick=2



; Look at errors
;resid_invert = REFORM(B_est1 - B)
resid_svd    = REFORM(reform(Bias_est) - Bias)
;plot, resid_invert, thick=2, title='Residual Error'
;oplot, resid_svd, color=4, thick=2


; Look at RMSE
;rmse_invert =  sqrt( TOTAL( reform(resid_invert)^2 ) / 280.0 )
rmse_svd    =  sqrt( TOTAL( reform(resid_svd)^2 ) / n )
print,  '  RMSE SVD =',rmse_svd

;window, 0
;plot, matts_wts
;window, 1
;plot, w_est2, color=2



; try the IDL regress function
;result = REGRESS(transpose(L), Bias)
;Bias_est_regress = transpose(L) ## reform(result)
;plot, bias_est_regress



END