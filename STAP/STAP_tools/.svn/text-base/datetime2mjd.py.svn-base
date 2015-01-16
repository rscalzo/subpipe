
import datetime

MJD_0=2400000.5
def datetime2mjd(dt):
    return datetime2jd(dt)-MJD_0

def mjd2datetime(mjd):
    return jd2datetime(mjd+MJD_0)

def datetime2jd(dt):
    yr=float(dt.year)
    mon=float(dt.month)
    day=float(dt.day)
    hr=float(dt.hour)
    min=float(dt.minute)
    sec=float(dt.second)
    msec=float(dt.microsecond)

    a = (14 - mon)//12
    y = yr+ 4800 - a
    m = mon + 12*a - 3
    j = day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045-0.5

    jd = j + (hr + (min + (sec/60.0))/60.0)/24.0+msec/1e6/3600./24.

    return jd

def jd2datetime(jd):
    j = int(jd+0.5)
    j = j+32044
    g = int(j/146097)
    dg = j %146097
    c = int((int(dg/36524) +1)*3/4)
    dc = dg - c *36524
    b = int(dc/1461)
    db = dc %1461
    a = int((int(db/365)+1)*3/4)
    da = db - a*365
    y = g*400+c*100+b*4+a
    m = int((da*5+308)/153)-2
    d = da - int((m+4)*153/5)+122

    yr = y-4800+int((m+2)/12)
    mon = (m+2)%12 +1
    day = d+1

    fracj = jd+0.5 -int(jd+0.5)
    hr = int(fracj*24)
    min = int((fracj*24-hr)*60)
    sec = int(((fracj*24-hr)*60-min)*60.)
    msec = int((((fracj*24-hr)*60-min)*60.-sec)*1e6)
    return datetime.datetime(yr,mon,day,hr,min,sec,msec)
