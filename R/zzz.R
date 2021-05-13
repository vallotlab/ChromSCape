.onLoad <- function(libname, pkgname) {
    set_ChromSCape_options()
    shiny::addResourcePath('www',
                        system.file('www',
                                    package = 'ChromSCape'))
}