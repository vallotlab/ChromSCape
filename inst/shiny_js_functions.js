/*!
 * JavaScript Shiny JS additional function
 * This script defines functions to save cookies / load cookies; able/disable tabs and all sort of function for improving user experience
 *
 *  See https://deanattali.com/shinyjs/extend and https://stackoverflow.com/questions/50049420/manipulate-existing-leaflet-map-in-a-shiny-app-with-javascript-using-shinyjs
 */

shinyjs.init_directory = function() {
                                          var my_cookie = Cookies.get("path"); 
                                          if(my_cookie){
                                             Shiny.onInputChange("path_cookie", my_cookie);/*Update value shiny variable path_cookie*/
                                          }
                                          console.log("shinyjs.init_directory");
                                        };

shinyjs.save_cookie = function(params){
                                            Cookies.set('path', params, { expires: 7 }); /*The cookie named 'path' will last 7 days*/
                                            var my_cookie = Cookies.get('path'); 
                                            Shiny.onInputChange("path_cookie", my_cookie);/*Update value shiny variable path_cookie*/
                                            console.log("shinyjs.save_cookie");
};  

shinyjs.disableTab = function(name) {
   /* var tab = $('.nav li a[data-value=' + name + ']');
  tab.bind('click.tab', function(e) {
    e.preventDefault();
    return false;
  });
  tab.addClass('disabled');
  console.log(tab); */
  prev = function(e) {
              e.preventDefault();
              return false;
            };
   var aNodeList = document.getElementsByTagName('a');
        for (var i = 0; i < aNodeList.length; i++) {
          if(aNodeList[i].getAttribute('data-value') == name) {
            
            $(aNodeList[i]).bind('click.tab', prev);
            $(aNodeList[i]).addClass('disabled');
          }
        }
        console.log("shinyjs.disableTab");
};

shinyjs.enableTab = function(name) {
  console.log("enabling");
     var aNodeList = document.getElementsByTagName('a');
        for (var i = 0; i < aNodeList.length; i++) {
          if(aNodeList[i].getAttribute('data-value') == name) {
            
            $(aNodeList[i]).unbind('click.tab');
            $(aNodeList[i]).removeClass('disabled');
          }
        }
        console.log("shinyjs.enableTab");
};
