echo ARCH = 'uname -m'

win = "i686"
mac = "darwin"

echo removing old eschaton
echo  

if [ "$(ARCH)" = 'win' ]; then
   rmdir -rf C:/python27/lib/site-packages/eschaton
fi

if [ "$(ARCH)" = 'mac' ]; then
   sudo rm -rf /Library/Python/2.7/site-packages/eschaton
fi

echo copying updated eschaton to site-packages'
echo  


if [ "$(ARCH)" = 'win' ]; then
   xcopy C:/python27/lib/site-packages/eschaton
fi

if [ "$(ARCH)" = 'mac' ]; then
   cp -R ~/GDrive/eschaton /Library/Python/2.7/site-packages/eschaton
fi

echo installing with setup.py

if [ "$(ARCH)" = 'win' ]; then
   python /Library/Python/2.7/site-packages/eschaton/setup.py install
fi

if [ "$(ARCH)" = 'mac' ]; then
   sudo python /Library/Python/2.7/site-packages/eschaton/setup.py install
fi