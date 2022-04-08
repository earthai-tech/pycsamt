

def poof_cli_usage (obj:str,
                    topmarker ='>',
                    bottommarker ='<',
                    mlen =77): 
    """
    Function focuses on display of the usage of each interpreter 
    command line (CLI). It uses the documentation of each command line
    module and give an illustrative example. It highlighted to instruct 
    the user how to introduce the CLI parameters 
    
    Parameters
    ----------
    obj : str 
        string module documentation.
    topmarker : str, optional
        Marker to introduce the section of usage. The default is '>'.
    bottommarker : str, optional
        Marker to end the section of usage. The default is '<'.
    mlen : int
        Marker length. The default is 77
    Returns
    -------
    string object to display
    
    Example
    -------
    >>> import pycsamt 
    >>> pycsamt.cli.poof_cli_usage (pycsamt.nm.__doc__)


    """ 
    usage_doc = [obj]
    usage_doc .insert(0, topmarker * mlen + '\n')
    usage_doc.append(bottommarker * (mlen + 5) )
    
    return ''.join(usage_doc)  