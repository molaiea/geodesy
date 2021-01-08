from django.shortcuts import render

# Create your views here.
def home(request) :
    return render(request, 'home.html')


def inverse(request) :
    return render(request, 'inverse.html')

def direct(request) :
    return render(request, 'direct.html')