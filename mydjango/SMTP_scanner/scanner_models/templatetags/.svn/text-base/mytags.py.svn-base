from django import template
import numpy as np

register = template.Library()
@register.inclusion_tag('table_header.html', takes_context=True)
def table_header(context, headers):
    return {
        'headers': headers,
        }

@register.filter
def get_tb_field(value,arg):
    return value[str(arg)]

@register.filter
def as_sixty(value,arg):
    if arg in "hour":
        usevalue=value/15.
        return '%02.0f %02.0f %04.1f'%sixty(usevalue)
    else:
        usevalue=value
        return '%02.0f %02.0f %02.0f'%sixty(usevalue)

def sixty(deg):
    adeg=abs(deg)
    dd=np.floor(adeg)
    mm=np.floor((adeg-dd)*60.)
    ss=((adeg-dd)*60.-mm)*60.
    if deg<0: dd= -1. *dd
    return (dd,mm,ss)
