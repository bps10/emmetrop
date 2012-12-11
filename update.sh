echo ARCH = 'ver'


echo removing old eschaton
echo.  

rmdir -rf C:/python27/lib/site-packages/eschaton


echo copying updated eschaton to site-packages
echo.

xcopy C:/python27/lib/site-packages/eschaton


echo installing with setup.py

python C:/python27/lib/site-packages/eschaton/setup.py install


