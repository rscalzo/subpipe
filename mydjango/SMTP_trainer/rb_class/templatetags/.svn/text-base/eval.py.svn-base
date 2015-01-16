from django import template
import re

register = template.Library()

@register.tag(name='eval')
def do_eval(parser, token):
    "Usage: {% eval %}1 + 1{% endeval %}"
    try:
        tag_name, arg = token.contents.split(None, 1)
    except ValueError:
        msg = '%r tag requires arguments' % token.contents[0]
        raise template.TemplateSyntaxError(msg)

    m = re.search(r'as (\w+)', arg)
    if m:
        var_name = m.groups()
    else:
        msg = '%r tag had invalid arguments' % tag_name
        raise template.TemplateSyntaxError(msg)

    nodelist = parser.parse(('endeval',))
    parser.delete_first_token()
    return EvalNode(nodelist,var_name)
    
class EvalNode(template.Node):
    def __init__(self,nodelist,var_name):
        self.nodelist=nodelist
        self.var_name=var_name
        
    def render(self, context):
        context[self.var_name]=eval(self.nodelist.render(context))
        return ''
