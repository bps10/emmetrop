echo 'removing old eschaton'
echo ' '

sudo rm -rf /Library/Python/2.7/site-packages/eschaton

echo 'copying updated eschaton to site-packages'
echo ' '

cp -R ~/GDrive/eschaton /Library/Python/2.7/site-packages/eschaton

echo 'installing with setup.py'

sudo python /Library/Python/2.7/site-packages/eschaton/setup.py install
