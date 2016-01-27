from django import template
import re
import datetime

register = template.Library()

@register.filter
def get_attr(value, arg):
    return getattr(value, arg)

@register.filter
def get_idx(value, idx):
    return value[idx]

@register.inclusion_tag('table_header.html', takes_context=True)
def table_header(context, headers):
    return { 'headers': headers, }

@register.filter
def get_tb_field(value,arg):
    return value[str(arg)]

@register.filter
def get_run(value):
    match=re.match('\w{3}(\d{8}_\d{6})',value)
    if match:
        return match.group(1)
    else:
        return None

@register.filter
def date2mjd(value,arg):
    date=value #datetime.datetime.strptime(value,'%Y-%m-%d')
    p0=date.toordinal()
    if isinstance(date,datetime.datetime): 
        p1=(date.hour+date.minute/60.+date.second/3600.)/24.
    else:
        p1=0.
    mjd=p0+p1-678575.5
    if arg==0:
        return '%d'%(int(mjd))
    else:
        return eval('"%.'+'%d'%int(arg)+'f"%(mjd)')

@register.filter
def mjd2date(value,arg):
    #arg controls how many digits after day
    #arg=0, integer day
    days=value+678575.5 #jd-1721425.0
    p0=int(days)
    p1=days-p0
    if arg==0:
        return datetime.datetime.fromordinal(p0).strftime('%Y-%m-%d')
    else:
        return datetime.datetime.fromordinal(p0).strftime('%Y-%m-%d')+eval('".%0'+'%d'%int(arg)+'d"%(round(p1*'+'%d'%(10**int(arg))+'))')


