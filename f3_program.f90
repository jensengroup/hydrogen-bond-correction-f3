program thirdgeneration_hbond_correction
    ! """
    ! Third Generation Hydrogen-bond Correction
    ! https://github.com/jensengroup/hydrogen-bond-correction-f3
    ! """
    implicit none

    character(len=*), parameter :: version = '1.0'

    ! Options
    character(len=32) :: arg
    logical :: param, debug, lgradient
    logical :: exists
    character(len=32) :: paramfile

    ! Model
    integer :: labels

    integer :: i


    param = .False.
    debug = .False.
    lgradient = .False.

    if(command_argument_count() == 0) then
        call print_help()
        stop
    endif

    do i = 1, command_argument_count()
        call get_command_argument(i, arg)

        select case (arg)
        case ('-g', '--gradient')
            lgradient = .True.

        case ('-d', '--debug')
            debug = .True.

        case ('-p', '--param')
            param = .True.

            call get_command_argument(i+1, paramfile)
            inquire(file=paramfile, exist=exists)
            if(.not.exists) then
                print '(a)', ''
                print '(2a)', 'Parameter file missing', paramfile
                print '(a)', ''
                stop
            endif

        case ('-v', '--version')
            print '(a)', ''
            print '(2a)', 'version ', version
            print '(a)', 'https://github.com/jensengroup/hydrogen-bond-correction-f3'
            call print_cite()
            stop

        case ('-h', '--help')
            call print_help()
            stop

        case default
            print '(a)', ''
            print '(2a)', 'Unrecognized option: ', arg
            print '(a)', ''
            stop

        end select
    enddo

    contains

    subroutine print_cite()
        print '(a)', ''
        print '(a)', 'Cite as:'
        print '(a)', 'Martin Korth'
        print '(a)', 'J. Chem. Theory Comput., 2010, 6 (12), pp 3808–3816'
        print '(a)', 'DOI: 10.1021/ct100408b'
        print '(a)', ''
    end subroutine print_cite

    subroutine print_help()
        print '(a)', ''
        print '(a)', 'usage: f3 [OPTIONS] <structure.xyz>'
        print '(a)', ''
        print '(a)', 'Calculate Hydrogen-bond correction for a XYZ structure'
        print '(a)', ''
        print '(a)', 'options:'
        print '(a)', ''
        print '(a)', '  -v, --version     print version information and exit'
        print '(a)', '  -h, --help        print usage information and exit'
        print '(a)', '  -p <file>, --param <file>   use parameter file <par.dat>'
        print '(a)', '  -d, --debug       print debug information'
        print '(a)', '  -g, --gradient    calculate and print gradient'
        print '(a)', ''
    end subroutine print_help


end program thirdgeneration_hbond_correction

