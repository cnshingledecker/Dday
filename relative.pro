; Change directory to the current location
; Note: The rest of the program assumes that you are in the sim. directory
CD, CURRENT=c & PRINT, c

; Put species to plot here
species = ["bH2"]
wrt     = ["bH2O"]
outfile = species[0] + "_rads.eps"


; Determine the number of species
nspecies = size(species,/N_ELEMENTS)

; Assign labels for plot legend
labels = ["Model","Experiment"]
;labels = species 
linestyles = INDGEN(nspecies)

; Plotting variables
LNTHCK = 3
charsize = 1.5


; Determine the number of file lines
abfile = c + '/ab/' + species[0] + ".ab"
nlines = FILE_LINES(abfile)

; Make the output array
; N+1 columns (N species + 1 time column)
; by
; nlines
outarray = FLTARR(nspecies,nlines-1)
wrtarray = FLTARR(nlines-1)


time = FLTARR(nlines-1)

; First get wrt species
abfile = c + '/ab/' + wrt[0] + ".ab"
temp = FLTARR(2)
OPENR, lun, abfile, /GET_LUN
; skip the first line
header = ''
READF, lun, header
FOR i = 0, nlines-2 DO BEGIN 
  READF, lun, temp,FORMAT='(E10.3,E13.6)'
  wrtarray[i] = temp[1]
ENDFOR
FREE_LUN, lun

; Now get the abundances of each species
FOREACH element, species, index DO BEGIN
  PRINT, "The species=",element," at index",index
  abfile = c + '/ab/' + species[index] + ".ab"
  temp = FLTARR(2)
  OPENR, lun, abfile, /GET_LUN
  ; skip the first line
  header = ''
  READF, lun, header
  FOR i = 0, nlines-2 DO BEGIN 
    READF, lun, temp,FORMAT='(E10.3,E13.6)'
    IF ( index EQ 0 ) THEN BEGIN
      time[i] = temp[0]
    ENDIF
    outarray[index,i] = temp[1]
  ENDFOR
  FREE_LUN, lun
ENDFOREACH

  maxtim = 0.0
  mintim = 0.0

  maxab = 0.0
  minab = 0.0

  xlab = "Time (s)"
;  xlab = "Fluence (10!E14!N ions/cm!E2!N)"
;  ylab = "[X] (10!E20!N cm!E-3!N)"
  ylab = "[" + species[0] + "]/[" + wrt[0] + "] %"

  FOR ii = 0,1 DO BEGIN
    IF (ii eq 0) THEN BEGIN
      SET_PLOT, 'PS'
      DEVICE, filename=outfile, /ENCAPSULATED
  ENDIF ELSE BEGIN
    DEVICE, /CLOSE
    SET_PLOT, 'X'
  ENDELSE

; Plot ranges ymin = min(outarray) ;1.0e-25 ; min(outarray)
ymin = min(outarray) - 0.5*min(outarray)
ymax = max(outarray) + 0.5*min(outarray)
xmin = min(time[*]) - (0.1-min(time[*]))
xmax = max(time[*])


print, ymax,ymin 



  ; Plot the ranges 
  PLOT,$
    time[*],$
    100.0*(outarray[0,*]/wrtarray[*]),$
    /XLOG,/YLOG,$
;    /YLOG,$
    linestyle=0,$ 
    xtitle=xlab,$ 
    ytitle=ylab,$
    yrange=[ymin,ymax],$
    xrange=[xmin,xmax],$
    XCHARSIZE=charsize,$
    YCHARSIZE=charsize,$
    XMARGIN=[13,5],$
    YMARGIN=[5.2,2], $
    XSTYLE=1,$
    YSTYLE=1,$
    /NODATA

  ; Plot the first species
  OPLOT,$
    time[*],$
    100*(outarray[0,*]/wrtarray[*]),$
    linestyle=0
    

  ; If there are more species, overplot these
  IF nspecies GT 1 THEN BEGIN
    FOR i = 1,nspecies-1 DO BEGIN
      PRINT, "Plotting ",species[i]
      OPLOT, time[*],100.0*(outarray[i,*]/wrtarray[*]),linestyle=i
    ENDFOR
  ENDIF


; Add a legend
PAT = BYTARR(10,10, /NOZERO)
PAT[5,5] = 255
LEGEND, labels, box=0, /LEFT_LEGEND;, pos=[0.72,0.35],/NORM,CHARSIZE=charsize

ENDFOR

input = 'a'
READ, input

IF input EQ '' THEN BEGIN
  WDELETE 
ENDIF ELSE BEGIN
  WDELETE 
ENDELSE


END



