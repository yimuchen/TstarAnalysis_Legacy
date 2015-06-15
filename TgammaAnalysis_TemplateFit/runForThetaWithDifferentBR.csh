#!/bin/csh

set templateName=TemplateFit_v21
set BRs="0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0"
#set BRs="0.5 0.6 0.7 0.8 0.9 1.0"
set BRs="0.8 1.0"

rm log_dummy temp_shape_validaton.C

foreach br($BRs)

echo $br
sed "s/float BR = 1.0/float BR = ${br}/g" ${templateName}.C | sed "s/${templateName}/temp_shape_validaton/g" >& temp_shape_validaton.C

root -l -b -q temp_shape_validaton.C++ >& log_dummy 

rm log_dummy temp_shape_validaton.C

end
