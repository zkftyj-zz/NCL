load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; This function is written to use bootstrap method to test the significance of correlation
; Y[*] is the parameter you want to test, time[*] is the time series, n[1] is the times of repeating calculating.
; Written by Kai Zhang

function bootstraptrend(X[*],Y[*],n[1])
local time,sloptmp,b0,Ytmp,bstorage,pcdf,bounds,findind,zpdf,k,genindice,ntime
begin
	ntime		= dimsizes(X)

	b0		= regline(X,Y)
	
	bstorage	= new((/n/),float)
	;generate bootstrap samples
	do	k = 0,n-1
		genindice	= toint(random_uniform(0,ntime-1,ntime))
		Ytmp		= Y(genindice)
		bstorage(k)= regline(X,Ytmp)
	end do

	opt		= True
	zpdf		= pdfx(bstorage,100,opt) ;100 bins

	pcdf		= cumsum(zpdf,0)
	bounds	= zpdf@bin_bounds
	findind	= ind(bounds.ge.b0)
	if (dimsizes(findind).gt.1) then
		return(pcdf(findind(0)))
	else
		return(pcdf(0))
	end if
end

