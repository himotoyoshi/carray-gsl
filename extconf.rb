require "mkmf"
require "carray/mkmf"

begin
  have_gsl = true
  $CPPFLAGS   << " " << `gsl-config --cflags`.chomp
  $LOCAL_LIBS << " " << `gsl-config --libs`.chomp
  `gsl-config --libs`.split(/\s+/).each do |opt|
    case opt
    when /\-l(.*)/
      unless have_library($1)
        have_gsl = false
        break
      end
    end
  end
rescue
end

$CFLAGS << " -Wall "  

if have_carray()
  if have_header("gsl/gsl_matrix.h") and have_gsl
    create_makefile("carray/carray_gsl")
  else
    open("Makefile", "w") { |io|
      io << "all:"       << "\n"
      io << "install:"   << "\n"        
      io << "clean:"     << "\n"
      io << "distclean:" << "\n"     
      io << "\trm -rf mkmf.log Makefile" << "\n"     
    }
  end
end
