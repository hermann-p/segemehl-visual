// Constructor for the BeamerLike function class
var BeamerLike = function () {
  this.SHOW_NAVBAR = true;
  this.SHOW_FOOTER = true;
  this.SHOW_HEADER = true;
  var navBar = document.createElement('header');
//  navBar.setAttribute('class', 'navBar');
//  navBar.setAttribute('id', 'navBar');
  navBar.setAttribute('style', 'visibility:hidden;position:fixed;top:0;height:2em;z-index:400;width:100%;text-align:center;')
  var insertionPoint = document.getElementsByClassName('slides');
  insertionPoint[0].parentNode.insertBefore(navBar, insertionPoint[0]);
  this.navBar = navBar;
}
// Apply configuration object

BeamerLike.prototype.configure = function (conf) {
  this.SHOW_NAVBAR = conf.show_navbar;
  this.SHOW_FOOTER = conf.show_footer;
  this.SHOW_HEADER = conf.show_header;
}
// Parse the document

BeamerLike.prototype.run = function () {
  var slideContainer = document.getElementsByClassName('slides')
  this.allSlides = slideContainer[0].getElementsByTagName('section')
  // loop through all slides after titleslide
  for (var ii = 1; ii < this.allSlides.length; ii++) {
    var currentSlide = this.allSlides[ii];
    // set data state for this slide
//    var oldAttr = currentSlide.getAttribute('data-state')+',' || '';
    currentSlide.setAttribute('data-state', 'header'+ii);
    // find content of slide's headline
    var heading = currentSlide.getElementsByTagName('h2');
    if (heading[0]) {
      var headline = heading[0].innerHTML;
      // define content for this data slide
      var headerInfo = document.createElement('style');
      headerInfo.innerHTML = '.header'+ii+' header:after{content:"' + headline + '";visibility:visible;}';
      currentSlide.insertBefore(headerInfo, currentSlide.firstChild);
      console.log(headline)
    }
  }
}
bl = new BeamerLike()
bl.run()