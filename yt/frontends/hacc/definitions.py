from numpy import ones_like

def fdark_matter(mask):
    mv = 2
    filt = ~((mask & (1<<mv)).astype('bool'))
    mv = 7
    filt = ~((mask & (1<<mv)).astype('bool'))
    return filt
    
def fbaryonic_gas(mask):
    mv = 2
    filt = (mask & (1<<mv)).astype('bool')
    mv = 3
    filt &= ~((mask & (1<<mv)).astype('bool'))
    mv = 4
    filt &= ~((mask & (1<<mv)).astype('bool'))
    mv = 5
    filt &= ~((mask & (1<<mv)).astype('bool'))
    mv = 7
    filt &= ~((mask & (1<<mv)).astype('bool'))
    mv = 8
    filt &= ~((mask & (1<<mv)).astype('bool'))
    return filt

def fstars(mask):
    mv = 4
    filt = (mask & (1<<mv)).astype('bool')
    mv = 7
    filt &= ~((mask & (1<<mv)).astype('bool'))
    return filt

def fagn(mask):
    mv = 3
    filt = (mask & (1<<mv)).astype('bool')
    mv = 7
    filt &= ~((mask & (1<<mv)).astype('bool'))
    mv = 8
    filt &= ~((mask & (1<<mv)).astype('bool'))
    mv = 9
    filt &= ~((mask & (1<<mv)).astype('bool'))
    return filt

def fwind(mask):
    mv = 5
    filt = (mask & (1<<mv)).astype('bool')
    mv = 7
    filt &= ~((mask & (1<<mv)).astype('bool'))
    return filt

def fsf_gas(mask):
    mv = 6
    filt = (mask & (1<<mv)).astype('bool')
    mv = 7
    filt &= ~((mask & (1<<mv)).astype('bool'))
    return filt
def fparticle(mask):
    return ones_like(mask).astype('bool')

FILTER_DEF = {
                # 'Particle':fparticle,
                'DarkMatter':fdark_matter, 
                'Gas':fbaryonic_gas,
                'Star':fstars,
                'AGN':fagn,
                'Wind':fwind,
                'SF_Gas':fsf_gas
                }
# maps simulation particle names to yt names
KEY_TO_YT = {
        'DarkMatter':'particle',
        'Gas':'gas',
        'Star':'star',
        'AGN':'agn',
        'Wind':'wind'
}
KEY_TO_SIM = {}
for k, v in KEY_TO_YT.items():
    KEY_TO_SIM[v] = k