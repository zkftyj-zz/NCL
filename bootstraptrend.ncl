load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; This function is written to use bootstrap method to test the significance of trend
; Y[*] is the parameter you want to test, time[*] is the time series, n[1] is the times of repeating calculating.
; Written by Kai Zhang

function bootstraptrend(Y[*],time[*],n[1])
local time,sloptmp,b0,Ytmp,bstorage,pcdf,bounds,findind,zpdf,k,genindice,i,j,ntime
begin
	ntime		= dimsizes(time)
	
	sloptmp	= new((/ntime-1,ntime-1/),float)

	Y1		= dtrend(Y,False)
	;b0 = median(xi-xj/i-j)
	do	i = 0,ntime-2
		do	j = 1,ntime-1
			sloptmp(i,j-1) = (Y1(j)-Y1(i))/(time(j)-time(i))
		end do
	end do
	b0		= dim_median(ndtooned(sloptmp))
	
	bstorage	= new((/n/),float)
	;generate bootstrap samples
	do	k = 0,n-1
		sloptmp	= 0
		genindice	= toint(random_uniform(0,ntime-1,ntime))
		Ytmp		= Y1(genindice)
		do	i = 0,ntime-2
			do	j = 1,ntime-1
			sloptmp(i,j-1) = (Ytmp(j)-Ytmp(i))/(time(j)-time(i))
			end do
		end do
		bstorage(k)= dim_median(ndtooned(sloptmp))
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

