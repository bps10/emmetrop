echo ARCH = 'uname -m'

mac = "darwin"

echo removing old eschaton
echo  

if [ "$(ARCH)" = 'mac' ]; then
   sudo rm -rf /Library/Python/2.7/site-packages/eschaton
fi

echo copying updated eschaton to site-packages
echo 


if [ "$(ARCH)" = 'mac' ]; then
   cp -R ~/GDrive/eschaton /Library/Python/2.7/site-packages/eschaton
fi

echo installing with setup.py


if [ "$(ARCH)" = 'mac' ]; then
   sudo python /Library/Python/2.7/site-packages/eschaton/setup.py install
fi