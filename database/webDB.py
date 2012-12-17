import urllib
import shlex
import os


baseURL = 'ftp://tofu.psych.upenn.edu/fulldb/'

def getUPENNimages(album, saveDir = None):
    if not saveDir:
        saveDir = album
    
    if os.path.isdir(saveDir):
        print 'found existing directory'
    else :
        os.makedirs(saveDir)
        print 'making new directory called {0}'.format(saveDir)
    
    print 'getting album names'
    print ' '    
    images = getAlbumNames(album, saveDir)
    
    print 'downloading in progress ... '
    print ' '
    for image in images:
        print image
        
        jpgURL = baseURL + album + '/' + image + '.JPG'
        urllib.urlretrieve(jpgURL, saveDir + '/' + image + '.JPG')
        
        auxURL = baseURL + album + '/' + image + '_AUX.mat'
        urllib.urlretrieve(auxURL, saveDir + '/' + image + '_AUX.mat')
    
    os.remove(saveDir + '/' + album + '.txt')
    
    print 'download to directory, {0}, complete'.format(saveDir)
    
def getAlbumNames(album, saveDir):
    albumURL = baseURL + album
    
    urllib.urlretrieve(albumURL, saveDir + '/' + album + '.txt')
    
    names = parseUpennFile(album, saveDir)
    
    return names

def parseUpennFile(album, saveDir):
    """
    
    """
    fil = open(saveDir + '/' + album + '.txt').read()
    foo = fil
    names = []
    partitioning = True

    while partitioning:
        
        f = foo.partition('\n')
        

        if f[0]:     
            parse = shlex.split(f[0])
            
            name = parse[8]
            
            if name[-3:].lower() == 'jpg':
                names.append(name[:-4])
            foo = f[2]
            partitioning = True
        else:
            partitioning = False
    
    return names        
    
if __name__ == '__main__':
    
    import argparse

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-g", "--group", type=str, default='cd01A',
                        help="indicate what directory to download")
    parser.add_argument("-d", "--dir", type=str, default='.',
                        help="indicate where to save files")
    
    args = parser.parse_args()
    getUPENNimages(args.group)
    