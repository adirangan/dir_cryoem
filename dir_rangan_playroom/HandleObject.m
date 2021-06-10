classdef HandleObject < handle
   properties
      O=[];
   end
 
   methods
      function obj=HandleObject(receivedObject)
         obj.O=receivedObject;
      end
   end
end
