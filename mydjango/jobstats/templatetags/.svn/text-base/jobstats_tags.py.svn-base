from django import template
import re

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
