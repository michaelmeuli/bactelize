cmake -DITK_DIR:PATH=/home/michael/Colocalization/ITK/build ..


first commited version worked with ITK master built on Apr 18 2014.
InsightToolkit-4.5.2 didn't work.



SCIFIO was fixed Thu Apr 10 12:45:05 2014 -0500: 

commit 9f6f245a10b7837e6d7303b872ed55611fd99b5a
Author: Mark Hiner <hinerm@gmail.com>
Date:   Thu Apr 10 12:45:05 2014 -0500

    ENH: bump to latest scifio-imageio
    
    Updated the SCIFIO-ImageIO to fix a pixel type detection error. It
    seemed that the returned pixel type was clashing with the
    UNKNOWNCOMPONENTTYPE constant, causing an unknown component error to be
    thrown erroneously.
    
    Change-Id: I8f691f82da9288709e48626dd13d102d9df15bff




