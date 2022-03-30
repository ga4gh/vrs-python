#Notes from Wes Goar on setup of vrs-python on a m1 mac
- @author Wes Goar
- date: 20220328
- macOS version: 12.3
- python version: 3.9.11

##Notes
- I installed seqrepo before trying this setup.
- I had already installed xcode and homebrew on my machine.

###Steps
1. brew install openssl
2. add the following statement in your .zshrc: export PATH="/opt/homebrew/opt/openssl@1.1/bin:$PATH"
3. add the following statements in your .zshenv:
4. Run brew info openssl to ascertain the correct environment path to put into the export statements

export LDFLAGS="-L/opt/homebrew/opt/openssl@1.1/lib" \
export CPPFLAGS="-I/opt/homebrew/opt/openssl@3/include" \
export PKG_CONFIG_PATH="/opt/homebrew/opt/openssl@1.1/lib/pkgconfig"

5. source ~/.zshrc
6. source ~/.zshenv
7. brew install libpq
8. add the following statement in your .zshrc: export PATH="/opt/homebrew/opt/libpq/bin:$PATH"
9. add the following statements in your .zshenv:
10. Run brew info libpq to ascertain the correct environment path to append to the appropriate flags:

export LDFLAGS="-L/opt/homebrew/opt/openssl@1.1/lib -L/opt/homebrew/opt/libpq/lib"
export CPPFLAGS="-I/opt/homebrew/opt/openssl@3/include -I/opt/homebrew/opt/libpq/include"
export PKG_CONFIG_PATH="/opt/homebrew/opt/openssl@1.1/lib/pkgconfig:/opt/homebrew/opt/libpq/lib/pkgconfig"

11. source ~/.zshenv
12. brew install postgres
13. add the following statement in your .zshrc: export PATH="/opt/homebrew/opt/postgresql@14/bin:$PATH"
    ###Make sure that you update the @14 with your own version
14. source ~/.zshrc
#Optional?
15. I also added this statement in my .zshrc but I am unsure if it was actually needed or not (Try without it first:
    export LIBRARY_PATH=$LIBRARY_PATH:/opt/homebrew/opt/openssl/lib/



